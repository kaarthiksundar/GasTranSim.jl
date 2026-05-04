

function explicit_staggered_grid_step!(ts::TransientSimulator, run_type::Symbol)
    advance_pipe_density_internal!(ts, run_type) # (n+1) level
    staggered_advance_junction_pressures!(ts, run_type) # pressure (n+1), flux (n+1/2)
    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_compute_pipe_end_fluxes_densities!, ts, key_array, run_type)
    advance_pipe_mass_flux_internal!(ts, run_type)
    return
end



function initialize_pipe_grid!(ts::TransientSimulator, ::Val{:explicit_staggered_grid})
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

function initialize_pipe_state!(ts::TransientSimulator, ::Val{:explicit_staggered_grid})
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
    end
    return
end

function advance_pipe_density_internal!(ts::TransientSimulator, run_type::Symbol)
    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_advance_pipe_density_internal!, ts, key_array, run_type)
    return
end

function advance_pipe_mass_flux_internal!(ts::TransientSimulator, run_type::Symbol)
    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_advance_pipe_mass_flux_internal!, ts, key_array, run_type)
    return
end



function staggered_advance_junction_pressures!(ts::TransientSimulator, _run_type::Symbol)
    x_node = get_density.(Ref(ts), form_nodal_pressure_vector(ts))
    residual_fun! = (r, x) -> assemble_junction_residual!(ts, x, r)
    Jacobian_fun! = (J, x) -> assemble_junction_Jacobian!(ts, x, J)

    x_node, converged, iter, res_norm = solve_newton_basic!(x_node, residual_fun!, Jacobian_fun!)
    
    converged || throw(DomainError(res_norm, "Newton solver did not converge for nodal densities"))

    check_limits(ts, x_node)

    # Update nodal pressures in ts.ref using the density solution x_node.
    for node_id = 1:length(x_node)
        p_val = get_pressure(ts, x_node[node_id])
        ref(ts, :node, node_id)["pressure_previous"] = ref(ts, :node, node_id)["pressure"]
        ref(ts, :node, node_id)["pressure"] = p_val
        ref(ts, :node, node_id)["is_updated"] = true
    end

    return
end

"""
    Advance the internal pipe densities using mass balance equation.
`` \\rho^{t+1}_i = \\rho^{t}_i + \\frac{\\Delta t}{\\Delta x} \\cdot \\left( \\phi^{t+}_{i} - \\phi^{t+}_{i+1} \\right)``
`` \\text{where, } t^+ = t + \\frac {\\Delta t} 2.``
"""
function _advance_pipe_density_internal!(ts::TransientSimulator, pipe_id::Int64)
    rho = ref(ts, :pipe, pipe_id)["density_profile"]
    phi = ref(ts, :pipe, pipe_id)["mass_flux_profile"]
    n = ref(ts, :pipe, pipe_id)["num_discretization_points"]
    dx = ref(ts, :pipe, pipe_id)["dx"]
    dt = params(ts, :dt)
    rho[2:(n-1)] = rho[2:(n-1)] + (dt / dx) * (phi[2:(n-1)] - phi[3:n])
    return
end

"""
    Function to compute \$\\phi^t_i\$ using densities
"""
function _invert_quadratic(a::Real, y::Real)::Real
    # can also write as 2 * y / (1 + sqrt( 1 + 4 * a * abs(y) ) )
    return sign(y) * (-1.0 + sqrt(1.0 + 4.0 * a * abs(y))) / (2.0 * a)
end

"""
    Advance pipe internal mass fluxes using momentum balance equation
``` a = \\frac{\\Delta t \\cdot \\beta}{\\rho^{t+1}_i + \\rho^{t+1}_{i+1}} ```
``` y = \\phi^t_i - \\frac{\\Delta t}{\\Delta x} \\cdot \\left( p^{t+1}_i - p^{t+1}_{i-1}\\right) - a \\cdot \\phi^t_i \\cdot |\\phi^t_i| ```
``` \\phi^{t+1}_i = \\operatorname{_invert_quadratic}(a, y) ```
"""
function _advance_pipe_mass_flux_internal!(ts::TransientSimulator, pipe_id::Int64)
    rho = ref(ts, :pipe, pipe_id)["density_profile"]
    phi = ref(ts, :pipe, pipe_id)["mass_flux_profile"]

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


    phi[2:n] = _invert_quadratic.(a_vec, y_vec)
    #update field
    ref(ts, :pipe, pipe_id)["fr_minus_mass_flux"] = phi[2]
    ref(ts, :pipe, pipe_id)["to_minus_mass_flux"] = phi[n]
    return
end

"""
    Compute fluxes and densities at both ends of the pipe
"""
function _compute_pipe_end_fluxes_densities!(ts::TransientSimulator, pipe_id::Int64)
    dx = ref(ts, :pipe, pipe_id)["dx"]  # can adjust with dx/2 for ghost point if needed
    dt = params(ts, :dt)
    n = ref(ts, :pipe, pipe_id)["num_discretization_points"]
    from_node_id = ref(ts, :pipe, pipe_id)["fr_node"]
    rho_prev = get_density(ts, ref(ts, :node, from_node_id)["pressure_previous"])
    rho = get_density(ts, ref(ts, :node, from_node_id)["pressure"])
    # at (n + 1/2) level
    ref(ts, :pipe, pipe_id)["fr_mass_flux"] =
        ref(ts, :pipe, pipe_id)["fr_minus_mass_flux"] + (rho - rho_prev) * (dx/dt)
    ref(ts, :pipe, pipe_id)["mass_flux_profile"][1] =
        ref(ts, :pipe, pipe_id)["fr_mass_flux"]
    # at (n + 1) level
    ref(ts, :pipe, pipe_id)["density_profile"][1] = rho
    to_node_id = ref(ts, :pipe, pipe_id)["to_node"]
    rho_prev = get_density(ts, ref(ts, :node, to_node_id)["pressure_previous"])
    rho = get_density(ts, ref(ts, :node, to_node_id)["pressure"])
    #at (n + 1/2) level
    ref(ts, :pipe, pipe_id)["to_mass_flux"] =
        ref(ts, :pipe, pipe_id)["to_minus_mass_flux"] + (rho_prev - rho) * (dx/dt)
    ref(ts, :pipe, pipe_id)["mass_flux_profile"][n+1] =
        ref(ts, :pipe, pipe_id)["to_mass_flux"]
    # at (n + 1) level
    ref(ts, :pipe, pipe_id)["density_profile"][n] = rho
    return
end



"""
    Assembles all pipe contributions to a node
"""
function _assemble_pipe_contributions_to_node(
    node_id::Int64,
    withdrawal::Real,
    ts::TransientSimulator)::Tuple{<:Real,<:Real}
    out_p = ref(ts, :outgoing_pipes, node_id)
    in_p = ref(ts, :incoming_pipes, node_id)
    jac_term = 0.0
    rhs_term = -1.0 * withdrawal # in input data, withdrawal is positive, but we want inflow positive
    pipe = ref(ts, :pipe)

    for p in out_p
        dx = pipe[p]["dx"] # can adjust with dx/2 for ghost point if needed
        jac_term += dx * pipe[p]["area"] / params(ts, :dt)
        rhs_term -= pipe[p]["fr_minus_mass_flux"] * pipe[p]["area"]
    end
    for p in in_p
        dx = pipe[p]["dx"] # can adjust with dx/2 for ghost point if needed
        jac_term += dx * pipe[p]["area"] / params(ts, :dt)
        rhs_term += pipe[p]["to_minus_mass_flux"] * pipe[p]["area"]
    end
    return jac_term, rhs_term
end



function assemble_junction_residual!(
    ts::TransientSimulator,
    x_node::AbstractArray,
    residual_node::AbstractArray,
)
    assemble_residual_for_nodes_WITHOUT_eqn_nos!(ts, x_node, residual_node)
    assemble_residual_for_nodes_WITH_eqn_nos!(ts, x_node, residual_node)
    return
end


function assemble_junction_Jacobian!(
    ts::TransientSimulator,
    x_node::AbstractArray,
    Jacobian::SparseMatrixCSC{Float64,Int64},
)
    assemble_Jacobian_for_nodes_WITHOUT_eqn_nos!(ts, x_node, Jacobian)
    assemble_Jacobian_for_nodes_WITH_eqn_nos!(ts, x_node, Jacobian)
    return
end


function assemble_residual_for_nodes_WITHOUT_eqn_nos!(ts::TransientSimulator, x_node::Vector{Float64}, residual_node::Vector{Float64})
    for (ci, compressor) in  get(ref(ts), :compressor, Dict())

        to_node = compressor["to_node"]
        from_node = compressor["fr_node"]
        ctrl_type, ctrl_val = control(ts, :compressor, ci, ref(ts, :current_time))
        if ctrl_type == c_ratio_control
            residual_node[ci] = get_pressure(ts, x_node[to_node]) - ctrl_val * get_pressure(ts, x_node[from_node]) 
        elseif ctrl_type == discharge_pressure_control
            residual_node[ci] = x_node[to_node] - get_density(ts, ctrl_val)
        elseif ctrl_type == flow_control
            # no contribution to Jacobian since flow is known
            continue
        end
    end
    return
end

# assemble compressor contribution to Jacobian 
function assemble_Jacobian_for_nodes_WITHOUT_eqn_nos!(ts::TransientSimulator, x_node::Vector{Float64}, J::SparseMatrixCSC{Float64,Int64})
    for (ci, compressor) in   get(ref(ts), :compressor, Dict())

        to_node = compressor["to_node"]
        from_node = compressor["fr_node"]
        ctrl_type, ctrl_val = control(ts, :compressor, ci, ref(ts, :current_time))
        if ctrl_type == c_ratio_control
            J[ci, from_node] = -ctrl_val * get_pressure_prime(ts, x_node[from_node])
            J[ci, to_node]  = get_pressure_prime(ts, x_node[to_node]) 
        elseif ctrl_type == discharge_pressure_control
            J[ci, to_node] = 1.0
        elseif ctrl_type == flow_control
            # no contribution to Jacobian since flow is known
            continue
        end
    end
    return
end


function assemble_residual_for_nodes_WITH_eqn_nos!(ts::TransientSimulator, x_node::Vector{Float64}, residual_node::Vector{Float64})
    for (node_id, node) in ref(ts, :node)

        eqn_num = ref(ts, :node, node_id)["eqn_number"]
        if  isnan(eqn_num)
            continue
        end
        ctrl_type, ctrl_val = control(ts, :node, node_id, ref(ts, :current_time))
        if ctrl_type == pressure_control
            residual_node[eqn_num] = x_node[node_id]- get_density(ts, ctrl_val)
            continue
        elseif ctrl_type == flow_control
            jac_term, rhs_pipe_term = _assemble_pipe_contributions_to_node(node_id, ctrl_val, ts)
            rhs_compressor_term = 
            _assemble_for_flow_control_compressors(node_id, ts)
            residual_node[eqn_num] += jac_term * x_node[node_id]  - jac_term * get_density(ts, ref(ts, :node, node_id, "pressure")) - rhs_pipe_term - rhs_compressor_term 
            continue
        end
    end
    return
end

function assemble_Jacobian_for_nodes_WITH_eqn_nos!(ts::TransientSimulator, x_node::Vector{Float64}, J::SparseMatrixCSC{Float64,Int64})
    for (node_id, node) in ref(ts, :node)

        eqn_num = ref(ts, :node, node_id)["eqn_number"]
        if  isnan(eqn_num)
            continue
        end
        ctrl_type, ctrl_val = control(ts, :node, node_id, ref(ts, :current_time))
        if ctrl_type == pressure_control
            J[eqn_num, node_id] = 1.0
            continue
        elseif ctrl_type == flow_control
            jac_term, _ = _assemble_pipe_contributions_to_node(node_id, ctrl_val, ts)
            J[eqn_num, node_id] += jac_term
            continue
        end
    end
    return
end

