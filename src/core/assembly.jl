

function assemble_compressor_equations!(ts::TransientSimulator, J::SparseMatrixCSC{Float64,Int64}, rhs::Vector{Float64})
    for (ci, compressor) in ref(ts, :compressor)
        to_node = compressor["to_node"]
        from_node = compressor["fr_node"]
        ctrl_type, ctrl_val = control(ts, :compressor, ci, ref(ts, :current_time))
        # right now doing rho_i - alpha * rho_j , should be doing
        # p(rho_i) - alpha * p(rho_j) and solve with NR
        if ctrl_type == c_ratio_control
            J[ci, from_node] = -1.0 * ctrl_val
            J[ci, to_node]  = 1.0
            rhs[ci] = 0.0
        elseif ctrl_type == discharge_pressure_control
            J[ci, to_node] = 1.0
            rhs[ci] = get_density(ts, ctrl_val)
        elseif ctrl_type == flow_control
            # no contribution to Jacobian since flow is known
            continue
        end
    end
    return
end
# assemble compressor contribution to Jacobian 

function assemble_nodal_equations!(ts::TransientSimulator, J::SparseMatrixCSC{Float64,Int64}, rhs::Vector{Float64})
    for (node_id, node) in ref(ts, :node)

        eqn_num = ref(ts, :node, node_id)["eqn_number"]
        # println(eqn_num)
        if  isnan(eqn_num)
            continue
        end
        ctrl_type, ctrl_val = control(ts, :node, node_id, ref(ts, :current_time))
        if ctrl_type == pressure_control
            J[eqn_num, node_id] = 1.0
            rhs[eqn_num] = get_density(ts, ctrl_val)
            continue
        elseif ctrl_type == flow_control
            jac_term, rhs_pipe_term = _assemble_pipe_contributions_to_node_new(node_id, ctrl_val, ts)
            rhs_compressor_term = _assemble_flow_control_compressor_contribution_to_node(node_id, ts)
            J[eqn_num, node_id] += jac_term
            rhs[eqn_num] += rhs_pipe_term + rhs_compressor_term + jac_term * get_density(ts, ref(ts, :node, node_id, "pressure"))
            continue
        end
    end
    return
end

"""
    Assembles all pipe contributions to a node
"""
function _assemble_pipe_contributions_to_node_new(
    node_id::Int64,
    withdrawal::Real,
    ts::TransientSimulator,
)::Tuple{<:Real,<:Real}
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

function _assemble_flow_control_compressor_contribution_to_node(
    node_id::Int64,
    ts::TransientSimulator,
)::Real
    out_c = ref(ts, :outgoing_compressors, node_id)
    in_c = ref(ts, :incoming_compressors, node_id)
    rhs_term = 0.0
    t = ref(ts, :current_time)
    for ci in in_c
        ctr, cmpr_val = control(ts, :compressor, ci, t)
        if ctr == flow_control
            rhs_term += cmpr_val # inflow positive
        end
    end
    for co in out_c
        ctr, cmpr_val = control(ts, :compressor, co, t)
        if ctr == flow_control
            rhs_term += (-1.0 * cmpr_val) # outflow negative
        end
    end
    return rhs_term
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
        (get_pressure(ts, rho[2:n]) - get_pressure(ts, rho[1:(n-1)])) -
        a_vec .* phi[2:n] .* abs.(phi[2:n])


    phi[2:n] = _invert_quadratic.(a_vec, y_vec)
    # update field
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


