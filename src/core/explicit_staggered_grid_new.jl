



function explicit_staggered_grid_step_new!(ts::TransientSimulator, run_type::Symbol)
    # precompute_block!(ts, advance_non_boundary_points!)
    advance_non_boundary_points!(ts)
    advance_junction_pressures!(ts, :explicit_staggered_grid_new, run_type)
    return
end

function initialize_pipe_grid!(ts::TransientSimulator, ::Val{:explicit_staggered_grid_new})
    for (key, pipe) in ref(ts, :pipe)
        # CFL condition c*dt/dx <= 0.9 => dx >= c*dt/0.9
        # with nondim dt, dx, we have nondim_dt/ nondim_dx < = 0.9 * mach_no
        c_inv = nominal_values(ts, :mach_num)
        num_segments =
            c_inv * (pipe["length"] * params(ts, :courant_number)) / params(ts, :base_dt)
        if num_segments < 1
            throw(CFLException(string(key)))
        end

        n = floor(Int64, num_segments) + 1
        ref(ts, :pipe, key)["num_discretization_points"] = n
        ref(ts, :pipe, key)["dx"] = pipe["length"] / (n - 1)
        ref(ts, :pipe, key)["density_profile"] = zeros(Float64, n)
        ref(ts, :pipe, key)["mass_flux_profile"] = zeros(Float64, n + 1)
    end
    return
end

function initialize_pipe_state!(ts::TransientSimulator, ::Val{:explicit_staggered_grid_new})
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
            pipe["fr_minus_mass_flux"] = initial_mass_flux # (dx/2)
            pipe["to_minus_mass_flux"] = initial_mass_flux # L-(dx/2)
            pipe["fr_mass_flux"] = initial_mass_flux # -(dx/2)
            pipe["to_mass_flux"] = initial_mass_flux # L+(dx/2)
        else
            flow_spl = initial_pipe_mass_flow(ts, key)
            pressure_spl = initial_pipe_pressure(ts, key)
            x_rho = LinRange(0, L, n)
            x_mid = x_rho[1:(n-1)] .+ dx / 2.0
            get_coeffs(flow_spl)[1]
            pipe["mass_flux_profile"] =
                [
                    get_coeffs(flow_spl)[1],
                    [flow_spl(x) for x in x_mid]...,
                    get_coeffs(flow_spl)[end],
                ] ./ area
            pipe["fr_minus_mass_flux"] = pipe["mass_flux_profile"][2]
            pipe["to_minus_mass_flux"] = pipe["mass_flux_profile"][end-1]
            pipe["fr_mass_flux"] = pipe["mass_flux_profile"][1]
            pipe["to_mass_flux"] = pipe["mass_flux_profile"][end]
            pipe["density_profile"][1:n] = [get_density(ts, pressure_spl(x)) for x in x_rho]
        end
        pipe["phi"] = zeros(n+1)
        pipe["rho"] = zeros(n)
    end
    return
end





"""
    Advance the internal pipe densities using mass balance equation.
`` \\rho^{t+1}_i = \\rho^{t}_i + \\frac{\\Delta t}{\\Delta x} \\cdot \\left( \\phi^{t+}_{i} - \\phi^{t+}_{i+1} \\right)``
`` \\text{where, } t^+ = t + \\frac {\\Delta t} 2.``
"""
function _advance_pipe_density_non_boundary_points!(ts::TransientSimulator, pipe_id::Int64)
    
    rho = ref(ts, :pipe, pipe_id)["rho"]
    phi = ref(ts, :pipe, pipe_id)["phi"]
    n = ref(ts, :pipe, pipe_id)["num_discretization_points"]
    dx = ref(ts, :pipe, pipe_id)["dx"]
    dt = params(ts, :dt)
    rho[2:(n-1)] = rho[2:(n-1)] + (dt / dx) * (phi[2:(n-1)] - phi[3:n]) # "rho" now updated
    return
end

"""
    Function to compute \$\\phi^t_i\$ using densities
"""
function _invert_quadratic_formula(a::Real, y::Real)::Real
    # can also write as 2 * y / (1 + sqrt( 1 + 4 * a * abs(y) ) )
    return sign(y) * (-1.0 + sqrt(1.0 + 4.0 * a * abs(y))) / (2.0 * a)
end

"""
    Advance pipe internal mass fluxes using momentum balance equation
``` a = \\frac{\\Delta t \\cdot \\beta}{\\rho^{t+1}_i + \\rho^{t+1}_{i+1}} ```
``` y = \\phi^t_i - \\frac{\\Delta t}{\\Delta x} \\cdot \\left( p^{t+1}_i - p^{t+1}_{i-1}\\right) - a \\cdot \\phi^t_i \\cdot |\\phi^t_i| ```
``` \\phi^{t+1}_i = \\operatorname{_invert_quadratic}(a, y) ```
"""
function _advance_pipe_mass_flux_non_boundary_points!(ts::TransientSimulator, pipe_id::Int64)
    rho = ref(ts, :pipe, pipe_id)["rho"]
    phi = ref(ts, :pipe, pipe_id)["phi"]

    n = ref(ts, :pipe, pipe_id)["num_discretization_points"]
    c = nominal_values(ts, :euler_num) / (nominal_values(ts, :mach_num))^2
    beta =
        ref(ts, :pipe, pipe_id, "friction_factor") /
        (2 * ref(ts, :pipe, pipe_id, "diameter"))
    a_vec = params(ts, :dt) * beta ./ (rho[2:n] + rho[1:(n-1)])
    y_vec =
        phi[2:n] -
        (c * params(ts, :dt) / ref(ts, :pipe, pipe_id, "dx")) *
        (get_pressure.(Ref(ts), rho[2:n]) - get_pressure.(Ref(ts), rho[1:(n-1)])) -
        a_vec .* phi[2:n] .* abs.(phi[2:n])


    phi[2:n] = _invert_quadratic_formula.(a_vec, y_vec) #"phi" now updated

    return
end

function _compute_pipe_end_fluxes!(ts::TransientSimulator, pipe_id::Int64, rho_from::T, rho_to::T) where {T<:Real}
    dx = ref(ts, :pipe, pipe_id)["dx"]  # can do dx/2 for ghost point if needed
    dt = params(ts, :dt)
    from_node_id = ref(ts, :pipe, pipe_id)["fr_node"]
    rho_prev_from = get_density(ts, ref(ts, :node, from_node_id)["pressure"])

    ref(ts, :pipe, pipe_id)["rho"][1] = rho_from
    ref(ts, :pipe, pipe_id)["rho"][end] = rho_to
    phi =  ref(ts, :pipe, pipe_id)["phi"]
    rho = ref(ts, :pipe, pipe_id)["rho"]


    phi[1] =
        phi[2] + (rho[1] - rho_prev_from) * (dx/dt)
    to_node_id = ref(ts, :pipe, pipe_id)["to_node"]
    rho_prev_to = get_density(ts, ref(ts, :node, to_node_id)["pressure"])
    #at (n + 1/2) level
    phi[end] =
        phi[end-1] + (rho_prev_to - rho[end]) * (dx/dt)
    return
end

function advance_non_boundary_points!(ts::TransientSimulator)

    for (pipe_id, pipe) in ref(ts, :pipe)
        ref(ts, :pipe, pipe_id)["rho"] = copy(ref(ts, :pipe, pipe_id)["density_profile"])
        ref(ts, :pipe, pipe_id)["phi"] = copy(ref(ts, :pipe, pipe_id)["mass_flux_profile"])
        _advance_pipe_density_non_boundary_points!(ts, pipe_id)
        _advance_pipe_mass_flux_non_boundary_points!(ts, pipe_id)
    end

    return
end
function _solve_pipe_state_staggered_grid!(ts::TransientSimulator,pipe_id::Int64,rho_from::T,rho_to::T,inertial_flag = zero(T))::Vector{T} where {T<:Real}

    area = T(ref(ts, :pipe, pipe_id)["area"])
    _compute_pipe_end_fluxes!(ts, pipe_id, rho_from, rho_to)
    phi = ref(ts, :pipe, pipe_id)["phi"]

    end_flows = T[area * phi[1], area * phi[end]]
    
    return end_flows

end
