function explicit_hyperbolic_step!(ts::TransientSimulator, run_type::Symbol)
    advance_junction_pressures!(ts, :explicit_hyperbolic, run_type)
    return
end


function initialize_pipe_grid!(ts::TransientSimulator, ::Val{:explicit_hyperbolic})
    dt_factor = 1.0
    while true
        try
            ts.params[:base_dt] = dt_factor * ts.params[:base_dt]
            for (key, pipe) in ref(ts, :pipe)
                # CFL condition c*dt/dx <= 0.9 => dx >= c*dt/0.9
                # with nondim dt, dx, we have nondim_dt/ nondim_dx < = 0.9 * mach_no
                 
                c_inv = nominal_values(ts, :mach_num)
                num_segments =
                    c_inv * (pipe["length"] * params(ts, :courant_number)) / params(ts, :base_dt)
                if num_segments < 2
                    throw(CFLException(string(key)))
                end
                n = floor(Int64, num_segments) + 1
                
                ref(ts, :pipe, key)["num_discretization_points"] = n
                ref(ts, :pipe, key)["dx"] = pipe["length"] / (n - 1)
                ref(ts, :pipe, key)["density_profile"] = zeros(Float64, n)
                ref(ts, :pipe, key)["mass_flux_profile"] = zeros(Float64, n)
            end
            break
        catch err

            if isa(err, CFLException)
                dt_factor *= 0.5
                println("CFL condition failed. Reducing Δt by half")
            else
                rethrow(err)
            end

        end
    end
    return
end

function initialize_pipe_state!(ts::TransientSimulator, ::Val{:explicit_hyperbolic})
    is_steady = false
    # We are assuming initial pipe pressures will not be provided only for steady
    # initial conditions. In the unsteady case we assume initial pipe pressures are given.
    if isempty(ts.initial_conditions[:pipe]["pressure"])
        @info "Pipes do not have initial pressure profile, will be computed assuming steady state flow"
        is_steady = true
    end
    for (key, pipe) in ref(ts, :pipe)
        area = pipe["area"]
        n = pipe["num_discretization_points"]
        dx = pipe["dx"]
        L = pipe["length"]
        fr_node = pipe["fr_node"]
        to_node = pipe["to_node"]
        if is_steady
            initial_mass_flux = initial_pipe_mass_flow(ts, key)(0.0) / area
            fill!(pipe["mass_flux_profile"], initial_mass_flux)
            initial_fr_pressure = ref(ts, :node, fr_node, "pressure")
            initial_to_pressure = ref(ts, :node, to_node, "pressure")
            density_at_first_sq = get_density(ts, initial_fr_pressure)^2
            density_at_last_sq = get_density(ts, initial_to_pressure)^2
            dL = dx / L
            pipe["density_profile"][1:n] = [
                sqrt(
                    density_at_last_sq * (i - 1) * dL + density_at_first_sq * (n - i) * dL,
                ) for i = 1:n
            ]
            pipe["fr_minus_mass_flux"] = initial_mass_flux # dx
            pipe["to_minus_mass_flux"] = initial_mass_flux # L-dx
            pipe["fr_mass_flux"] = initial_mass_flux # 0
            pipe["to_mass_flux"] = initial_mass_flux # L
        else
            flow_spl = initial_pipe_mass_flow(ts, key)
            pressure_spl = initial_pipe_pressure(ts, key)
            x_arr = LinRange(0, L, n)
            get_coeffs(flow_spl)[1]
            pipe["mass_flux_profile"][1:n] =
                [flow_spl(x) for x in x_arr] ./ area
            pipe["fr_minus_mass_flux"] = pipe["mass_flux_profile"][2]
            pipe["to_minus_mass_flux"] = pipe["mass_flux_profile"][n-1]
            pipe["fr_mass_flux"] = pipe["mass_flux_profile"][1]
            pipe["to_mass_flux"] = pipe["mass_flux_profile"][end]
            pipe["density_profile"][1:n] = [get_density(ts, pressure_spl(x)) for x in x_arr]
        end
        pipe["phi"] = zeros(n)
        pipe["rho"] = zeros(n)
    end
    return
end

function _maccormack_step!(ts::TransientSimulator,
    pipe_id::Int64,
    rho_from::T,
    rho_to::T;
    inertial_flag::T = zero(T)) where {T<:Real}

    n = ref(ts, :pipe, pipe_id)["num_discretization_points"]
    @assert n > 2 "McCormack step requires at least 3 spatial points"

    rho = T.(ref(ts, :pipe, pipe_id)["density_profile"]) #T. makes it a copy
    phi = T.(ref(ts, :pipe, pipe_id)["mass_flux_profile"])
    dx = T(ref(ts, :pipe, pipe_id)["dx"])
    dt = T(params(ts, :dt))
    nondim = T(nominal_values(ts, :euler_num) / (nominal_values(ts, :mach_num))^2)
    beta = T(
        ref(ts, :pipe, pipe_id, "friction_factor") /
        (2 * ref(ts, :pipe, pipe_id, "diameter")),
    )
    grav_coeff = T(
        ref(ts, :pipe, pipe_id, "sin_incline") / (nominal_values(ts, :froude_num))^2,
    )

    # Sources at current state.
    S_rho = zeros(T, n)
    S_phi = -beta .* phi .* abs.(phi) ./ rho .+ grav_coeff .* rho

    # Fluxes at current state.
    F_rho = copy(phi)
    F_phi = inertial_flag .* (phi .^ 2 ./ rho) .+ nondim .* get_pressure.(Ref(ts), rho)

    # Predictor.
    rho_p = copy(rho)
    phi_p = copy(phi)
    rho_p[1:n-1] .= rho[1:n-1] .- (dt / dx) .* (F_rho[2:n] .- F_rho[1:n-1]) .+ dt .* S_rho[1:n-1]
    phi_p[1:n-1] .= phi[1:n-1] .- (dt / dx) .* (F_phi[2:n] .- F_phi[1:n-1]) .+ dt .* S_phi[1:n-1]

    # Predicted-state sources and fluxes.
    S_rho_p = zeros(T, n)
    S_phi_p = -beta .* phi_p .* abs.(phi_p) ./ rho_p .+ grav_coeff .* rho_p
    F_rho_p = copy(phi_p)
    F_phi_p = inertial_flag .* (phi_p .^ 2 ./ rho_p) .+ nondim .* get_pressure.(Ref(ts), rho_p)

    # Corrector (backward difference).
    rho[2:n-1] .= T(0.5) .* (
        rho[2:n-1] .+ rho_p[2:n-1] .-
        (dt / dx) .* (F_rho_p[2:n-1] .- F_rho_p[1:n-2]) .+
        dt .* S_rho_p[2:n-1]
    )
    phi[2:n-1] .= T(0.5) .* (
        phi[2:n-1] .+ phi_p[2:n-1] .-
        (dt / dx) .* (F_phi_p[2:n-1] .- F_phi_p[1:n-2]) .+
        dt .* S_phi_p[2:n-1]
    )

    # Network-coupled boundary conditions: nodal densities on both ends.
    rho[1] = rho_from
    rho[n] = rho_to

    # Extrapolate edge fluxes.
    phi[1] = 2 * phi[2] - phi[3]
    phi[n] = 2 * phi[n-1] - phi[n-2]

    if T == Float64
        ref(ts, :pipe, pipe_id)["rho"] = rho
        ref(ts, :pipe, pipe_id)["phi"] = phi
    end

    return rho, phi
end

function _solve_pipe_state_maccormack!(ts::TransientSimulator,pipe_id::Int64,rho_from::T,rho_to::T,inertial_flag = zero(T))::Vector{T} where {T<:Real}

    area = T(ref(ts, :pipe, pipe_id)["area"])

    _, phi = _maccormack_step!(ts, pipe_id, rho_from, rho_to; inertial_flag = inertial_flag)

    end_flows = T[area * phi[1], area * phi[end]]
    
    return end_flows

end


