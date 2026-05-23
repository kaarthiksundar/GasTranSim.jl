function implicit_parabolic_step!(ts::TransientSimulator, run_type::Symbol)
    advance_junction_pressures!(ts, :implicit_parabolic, run_type) # pressure (n+1), flux (n+1/2)
    # update density and mass_flux profile
    # update end fluxes
    return
end


function initialize_pipe_grid!(ts::TransientSimulator, ::Val{:implicit_parabolic})
    for (key, pipe) in ref(ts, :pipe)
        pipe_segments = get(params(ts), :pipe_segments, 20)
        n = pipe_segments + 1
        ref(ts, :pipe, key)["num_discretization_points"] = n
        ref(ts, :pipe, key)["dx"] = pipe["length"] / (n - 1)
        ref(ts, :pipe, key)["density_profile"] = zeros(Float64, n)
        ref(ts, :pipe, key)["mass_flux_profile"] = zeros(Float64, n)
    end
    return
end

function initialize_pipe_state!(ts::TransientSimulator, ::Val{:implicit_parabolic})
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

function _solve_pipe_state_parabolic!(
    ts::TransientSimulator,
    pipe_id::Int64,
    rho_from::T,
    rho_to::T)::Vector{T} where {T<:Real}
    n = ref(ts, :pipe, pipe_id)["num_discretization_points"]
    x = zeros(T, 2 * n)
    x[1:n] = T.(ref(ts, :pipe, pipe_id)["density_profile"])
    x[n+1:2*n] = T.(ref(ts, :pipe, pipe_id)["mass_flux_profile"])
    rho_old = copy(x[1:n])
    area = T(ref(ts, :pipe, pipe_id)["area"])

    residual_fun! = (r, x) -> _pipe_residual!(
        r,
        ts,
        pipe_id,
        x,
        rho_old,
        rho_from,
        rho_to
    )
    Jacobian_fun! = (J, x) -> _pipe_jacobian!(
        J,
        ts,
        pipe_id,
        x
    )


    x, converged, _, res_norm = solve_newton_basic!(x, residual_fun!, Jacobian_fun!)

    
    converged || throw(DomainError(res_norm, "Newton solver did not converge for pipe $pipe_id"))

    if T == Float64
        ref(ts, :pipe, pipe_id)["rho"] = x[1:n]
        ref(ts, :pipe, pipe_id)["phi"] = x[n+1:2*n]
        ref(ts, :pipe, pipe_id)["fr_mass_flux"] = x[n+1]
        ref(ts, :pipe, pipe_id)["to_mass_flux"] = x[2*n]
    end

    end_flows = T[area * x[n+1], area * x[2*n]]
    
    return end_flows

end

"""
Residual for implicit parabolic pipe solve, matching solve_pipe_state! layout.

Unknown/state ordering:
- x[1:n]     -> rho
- x[n+1:2n]  -> phi

Residual ordering:
- r[1:n]     -> mass equations
- r[n+1:2n]  -> momentum equations
"""
function _pipe_residual!(
    r::Vector{T},
    ts::TransientSimulator,
    pipe_id::Int64,
    x::Vector{T},
    rho_old::Vector{T},
    rho_from::T,
    rho_to::T) where {T<:Real}

    n = div(length(x) , 2)

    rho = x[1:n]
    phi = x[n+1:2*n]
    @assert length(phi) == n
    @assert length(rho_old) == n
    @assert length(r) == 2n
    @assert n > 2

    fill!(r, zero(T))

    dx = ref(ts, :pipe, pipe_id)["dx"]
    dt = params(ts, :dt)

    # -------------------------
    # Mass balance residual
    # -------------------------
    r[1] = rho[1] - rho_from
    r[2:n-1] .= rho[2:n-1] .- rho_old[2:n-1] .+ (dt / dx) .* (phi[2:n-1] .- phi[1:n-2])
    r[n] = rho[n] - rho_to

    # -------------------------
    # Momentum balance residual
    # -------------------------
    rho_x = zeros(T, n)
    rho_x[1] = (rho[2] - rho[1]) / dx
    rho_x[2:n-1] .= (rho[3:n] .- rho[1:n-2]) ./ (2.0 * dx)
    rho_x[n] = (rho[n] - rho[n-1]) / dx

    phi_sq = phi .^ 2
    conv_term = zeros(T, n)
    conv_term[1] = (phi_sq[2] / rho[2] - phi_sq[1] / rho[1]) / dx
    conv_term[2:n-1] .= (phi_sq[3:n] ./ rho[3:n] .- phi_sq[1:n-2] ./ rho[1:n-2]) ./ (2.0 * dx)
    conv_term[n] = (phi_sq[n] / rho[n] - phi_sq[n-1] / rho[n-1]) / dx

    c2 = get_pressure_prime.(Ref(ts), rho)  # vectorized in your types.jl
    beta = ref(ts, :pipe, pipe_id, "friction_factor") /
        (2 * ref(ts, :pipe, pipe_id, "diameter"))
    nondim = nominal_values(ts, :euler_num) / (nominal_values(ts, :mach_num))^2
    phi_fun = phi .* abs.(phi)
    inertial_switch = params(ts, :inertial_flag) ? 1 : 0
    

    r[n+1:2n] .=
        nondim .* c2 .* rho_x .+
        inertial_switch .* conv_term .+
        beta .* phi_fun ./ rho .-
        ref(ts, :pipe, pipe_id, "sin_incline") .* rho / (nominal_values(ts, :froude_num))^2

    return r
end


"""
Build Jacobian J = d r / d x for x = [rho; phi], where
r[1:n]     = mass-balance residual
r[n+1:2n]  = momentum-balance residual
matching solve_pipe_state! residual definitions.
"""
function _pipe_jacobian!(
    J::AbstractMatrix{T},
    ts::TransientSimulator,
    pipe_id::Int64,
    x::Vector{T}) where {T<:Real}

    n = div(length(x), 2)
    rho = x[1:n]
    phi = x[n+1:2*n]
    @assert length(phi) == n
    @assert size(J, 1) == 2n && size(J, 2) == 2n

    if J isa SparseMatrixCSC
        fill!(J.nzval, zero(T))
    else
        fill!(J, zero(T))
    end

    dx = ref(ts, :pipe, pipe_id)["dx"]

    rho_col(i) = i
    phi_col(i) = n + i
    mom_row(i) = n + i

    # -------------------------
    # Mass rows: r[1:n]
    # -------------------------
    J[1, rho_col(1)] = 1.0
    for i in 2:n-1
        J[i, rho_col(i)] = 1.0
        J[i, phi_col(i)] = 1.0 / dx * params(ts, :dt)
        J[i, phi_col(i - 1)] = -1.0 / dx * params(ts, :dt)
    end
    J[n, rho_col(n)] = 1.0

    # c2 = dp/drho, c2p = d2p/drho2
    c2 = get_pressure_prime.(Ref(ts), rho)
    c2p = get_pressure_double_prime.(Ref(ts), rho)
    nondim = nominal_values(ts, :euler_num) / (nominal_values(ts, :mach_num))^2
    inertial_switch = params(ts, :inertial_flag) ? 1 : 0
    


    # helper for conv-term contributions from g_k = phi_k^2 / rho_k
    function add_conv_g!(row::Int, k::Int, coeff::Real)
        J[row, rho_col(k)] += inertial_switch * coeff * (-(phi[k]^2) / (rho[k]^2))
        J[row, phi_col(k)] += inertial_switch * coeff * (2.0 * phi[k] / rho[k])
    end

    for i in 1:n
        row = mom_row(i)

        # rho_x(i) and its rho-derivatives
        rho_x_i = zero(T)
        if i == 1
            rho_x_i = (rho[2] - rho[1]) / dx
            J[row, rho_col(1)] += nondim * c2[i] * (-1.0 / dx)
            J[row, rho_col(2)] += nondim * c2[i] * ( 1.0 / dx)
        elseif i == n
            rho_x_i = (rho[n] - rho[n - 1]) / dx
            J[row, rho_col(n - 1)] += nondim * c2[i] * (-1.0 / dx)
            J[row, rho_col(n)]     += nondim * c2[i] * ( 1.0 / dx)
        else
            rho_x_i = (rho[i + 1] - rho[i - 1]) / (2.0 * dx)
            J[row, rho_col(i - 1)] += nondim * c2[i] * (-1.0 / (2.0 * dx))
            J[row, rho_col(i + 1)] += nondim * c2[i] * ( 1.0 / (2.0 * dx))
        end

        # chain term from c2(rho_i) * rho_x(i)
        J[row, rho_col(i)] += nondim * c2p[i] * rho_x_i

        # inertial conv_term
        if i == 1
            add_conv_g!(row, 2,  1.0 / dx)
            add_conv_g!(row, 1, -1.0 / dx)
        elseif i == n
            add_conv_g!(row, n,      1.0 / dx)
            add_conv_g!(row, n - 1, -1.0 / dx)
        else
            add_conv_g!(row, i + 1,  1.0 / (2.0 * dx))
            add_conv_g!(row, i - 1, -1.0 / (2.0 * dx))
        end

        # friction term: beta * phi_i * abs(phi_i) / rho_i
        beta = ref(ts, :pipe, pipe_id, "friction_factor") /
        (2 * ref(ts, :pipe, pipe_id, "diameter"))
        J[row, rho_col(i)] += -beta * phi[i] * abs(phi[i]) / (rho[i]^2)
        J[row, phi_col(i)] +=  beta * (2.0 * abs(phi[i])) / rho[i]  # subgradient at 0 -> 0

        # gravity term: -G * rho_i
        J[row, rho_col(i)] += -ref(ts, :pipe, pipe_id, "sin_incline") / (nominal_values(ts, :froude_num))^2
    end

    return J
end


function implicit_advance_junction_pressures!(ts::TransientSimulator, _run_type::Symbol)
    x_node = get_density.(Ref(ts), form_nodal_pressure_vector(ts))
    
    problem_fun! = (r, J, x) -> assemble_network_problem!(ts, x, r, J)

    x_node, converged, iter, res_norm = NR_solve!(x_node, problem_fun!)
    # println(iter)
    converged || throw(DomainError(res_norm, "Newton solver did not converge for nodal densities"))

    check_limits(ts, x_node)

    # Update nodal pressures in ts.ref using the density solution x_node.
    for node_id = 1:length(x_node)
        p_val = get_pressure(ts, x_node[node_id])
        ref(ts, :node, node_id)["pressure_previous"] = ref(ts, :node, node_id)["pressure"]
        ref(ts, :node, node_id)["pressure"] = p_val
        ref(ts, :node, node_id)["is_updated"] = true
    end

    # update mass flux and density profiles in pipes
    for (pipe_id, pipe) in ref(ts, :pipe)
        pipe["density_profile"] = pipe["rho"]
        pipe["mass_flux_profile"] = pipe["phi"]
        pipe["fr_mass_flux"] = pipe["phi"][1]
        pipe["to_mass_flux"] = pipe["phi"][end]
    end

    return
end

