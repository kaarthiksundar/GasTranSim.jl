
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
    rho[2:n-1] = rho[2:n-1] + (dt / dx) *  (phi[2:n-1] - phi[3:n])
    return
end


function _assemble_eqns!(ts::TransientSimulator, pipe_id::Int64)
    n = ts.ref[:pipe][pipe_id]["num_discretization_points"]
    rhs = zeros(Float64, 2*n) 

    i_vec, j_vec, k_vec =  Vector{Int32}(), Vector{Int32}(), Vector{Float64}()
    sizehint!(i_vec, 10*n)
    sizehint!(j_vec, 10*n)
    sizehint!(k_vec, 10*n)
    for ele = 1:n-1  # elements
        # 4 eqns 
        push!(i_vec, ele)  
        push!(j_vec, 2*ele-1)
        push!(k_vec, 1)

        push!(i_vec, ele) 
        push!(j_vec, 2*ele+1)
        push!(k_vec, 1)

        push!(i_vec, ele)  
        push!(j_vec, 2*ele+2)
        push!(k_vec, 1)

        push!(i_vec, ele)  # 4 eqns 
        push!(j_vec, 2*ele)
        push!(k_vec, 1)

        rhs[ele] = ts.ref[:pipe][pipe_id]["pressure_profile"][2*ele-1] + ts.ref[:pipe][pipe_id]["pressure_profile"][2*ele+1] + 
    end

    

"""
    Sets the pressure and density at a given node to a value
"""
function _set_pressure_at_node!(node_id::Int64, pressure::Real, ts::TransientSimulator)
    (ref(ts, :node, node_id)["is_updated"] == true) && (return)

    ref(ts, :node, node_id)["pressure_previous"] = ref(ts, :node, node_id)["pressure"]
    ref(ts, :node, node_id)["density_previous"] = ref(ts, :node, node_id)["density"]

    ref(ts, :node, node_id)["pressure"] = pressure
    ref(ts, :node, node_id)["density"] = get_density(ts, pressure)

    ref(ts, :node, node_id)["is_updated"] = true
    return
end

"""
    Sets pressure and density at node across a compressor
"""
function _set_pressure_at_node_across_compressors!(node_id::Int64, pressure::Real, ts::TransientSimulator)
    t = ref(ts, :current_time)
    for ci in ref(ts, :incoming_compressors, node_id)
        ctrl_type, val = control(ts, :compressor, ci, t)
        i = ref(ts, :compressor, ci, "fr_node")
        if ctrl_type == c_ratio_control
            _set_pressure_at_node!(i, pressure/val, ts)
        elseif ctrl_type == discharge_pressure_control
            continue # pressure at both ends must be known by now
        elseif ctrl_type == flow_control
            continue
        end
    end
    for ci in ref(ts,:outgoing_compressors, node_id)
        ctrl_type, val = control(ts, :compressor, ci, t)
        i = ref(ts, :compressor, ci, "to_node")
        if ctrl_type == c_ratio_control
            _set_pressure_at_node!(i, pressure * val, ts)
        elseif ctrl_type == discharge_pressure_control
            _set_pressure_at_node!(i, val, ts)
        elseif ctrl_type == flow_control
            continue
        end
    end
    return
end

"""
    Solves for the pressure at the given node (with flow control) and its neighbors
"""
function _solve_for_pressure_at_node_and_neighbours!(node_id::Int64, withdrawal::Real, ts::TransientSimulator)
    t = ref(ts, :current_time)
    in_c = ref(ts, :incoming_compressors, node_id)
    index = findfirst(ci -> control(ts, :compressor, ci, t)[1] == discharge_pressure_control, in_c)

    # we know there can be at most one such compressor
    if isa(index, Nothing)
        # _calculate_pressure_for_node_without_incoming_discharge_pressure_control!(node_id, withdrawal, ts)
        out_p = ref(ts, :outgoing_pipes, node_id)
        in_p = ref(ts, :incoming_pipes, node_id)
        in_c = ref(ts, :incoming_compressors, node_id)
        out_c = ref(ts, :outgoing_compressors, node_id)
        vertex_vol = 0.0
        term1 = 0.0
        term2 = -1.0 * withdrawal # in input data, withdrawal is positive, but we want inflow positive
        for p in out_p
            pipe = ref(ts, :pipe)
            term1 += -1 * pipe[p]["fr_mass_flux"] * pipe[p]["area"]
            vertex_vol += pipe[p]["area"] * ts.params[:dx]
        end
        for p in in_p
            pipe = ref(ts, :pipe)
            term1 += pipe[p]["to_mass_flux"] * pipe[p]["area"]
            vertex_vol += pipe[p]["area"] * ts.params[:dx]
        end
        for co in out_c
            comp = ref(ts, :compressor)
            term1 += -1 * comp[co]["flow"] 
        end
        for ci in in_c
            comp = ref(ts, :compressor)
            term1 += comp[ci]["flow"] 
        end
        rho_new  =  ref(ts, :node, node_id)["density_previous"] + ts.params[:dt] * (term1 + term2) / vertex_vol
        p_new = get_pressure(ts, rho_new)
        _set_pressure_at_node!(node_id, p_new, ts)
        _set_pressure_at_node_across_compressors!(node_id, p_new, ts)


    else
        _calculate_pressure_for_node_with_incoming_discharge_pressure_control!(in_c[index], ts)
       
    end
    return
end

"""
    Calculate and update the pressures at nodes sans incoming discharge pressure-controlled compressor
"""
# function _calculate_pressure_for_node_without_incoming_discharge_pressure_control!(node_id::Int64, withdrawal::Real, ts::TransientSimulator)
#     s  = _assemble_compressor_contributions_to_node_without_incoming_discharge_pressure_control!(node_id, ts)
#     if isa(s, Nothing)
#         # means vertex pressure is already set
#         pressure = ref(ts, :node, node_id, "pressure")
#         _set_pressure_at_node_across_compressors!(node_id, pressure, ts)
#         return
#     end
#     s1, s2 = s
#     t1, t2 = _assemble_pipe_contributions_to_node(node_id, withdrawal, 1.0, ts)
#     rho_prev = get_density(ts, ref(ts, :node, node_id, "pressure"))
#     rho = rho_prev + (t2 + s2) / (t1 + s1)
#     pressure = get_pressure(ts, rho)
#     _set_pressure_at_node!(node_id, pressure, ts)
#     _set_pressure_at_node_across_compressors!(node_id, pressure, ts)
#     return
# end


"""
    Calculate and update the pressure at vertex with single compressor with flow_control 
"""
# function _set_pressure_for_node_with_single_flow_control_compressor!(node_id::Int64, net_withdrawal::Real, ts::TransientSimulator)
#     t1, t2 = _assemble_pipe_contributions_to_node(node_id, -net_withdrawal, 1.0, ts)
#     rho_prev = get_density(ts, ref(ts, :node, node_id, "pressure"))
#     rho = rho_prev + (t2  / t1)
#     pressure = get_pressure(ts, rho)
#     _set_pressure_at_node!(node_id, pressure, ts)
#     return
# end


"""
    Calculate and update the pressures at nodes with incoming discharge pressure controlled compressor
"""
function _calculate_pressure_for_node_with_incoming_discharge_pressure_control!(compressor_id::Int64, ts::TransientSimulator)
    t = ref(ts, :current_time)
    # node_id = ref(ts, :compressor, compressor_id)["fr_node"]
    base_node_id = ref(ts, :compressor, compressor_id)["to_node"]
    ctrl, discharge_pressure = control(ts, :compressor, compressor_id, t)
    @assert ctrl == discharge_pressure_control
    _set_pressure_at_node!(base_node_id, discharge_pressure, ts)
    _set_pressure_at_node_across_compressors!(base_node_id, discharge_pressure, ts)
    return
end

"""
    Assembles all pipe contributions to a node
"""
function _assemble_pipe_contributions_to_node(node_id::Int64, withdrawal::Real, mult_factor::Float64, ts::TransientSimulator)::Tuple{<:Real, <:Real}
    out_p = ref(ts, :outgoing_pipes, node_id)
    in_p = ref(ts, :incoming_pipes, node_id)
    term1 = 0.0
    term2 = -1.0 * withdrawal # in input data, withdrawal is positive, but we want inflow positive
    pipe = ref(ts, :pipe)
    for p in out_p
        term1 += mult_factor * pipe[p]["dx"] * pipe[p]["area"] / params(ts, :dt)
        term2 -= pipe[p]["fr_minus_mass_flux"] * pipe[p]["area"]
    end
    for p in in_p
        term1 += mult_factor * pipe[p]["dx"] * pipe[p]["area"] / params(ts, :dt)
        term2 += pipe[p]["to_minus_mass_flux"] * pipe[p]["area"]
    end
    return term1, term2
end

"""
    Assembles all comppressor contributions to a node sans incoming discharge pressure-controlled compressors
"""
function _assemble_compressor_contributions_to_node_without_incoming_discharge_pressure_control!(node_id::Int64, ts::TransientSimulator)::Union{Nothing, Tuple{<:Real, <:Real}}
    out_c = ref(ts, :outgoing_compressors, node_id)
    in_c = ref(ts, :incoming_compressors, node_id)
    term1 = 0.0
    term2 = 0.0
    p_prev = ref(ts, :node, node_id, "pressure")
    rho_prev = get_density(ts, p_prev)
    t = ref(ts, :current_time)
    for ci in in_c
        ctr, cmpr_val = control(ts, :compressor, ci, t)
        node_across_ci_id = ref(ts, :compressor, ci)["fr_node"]
        ctrl_type, val = control(ts, :node, node_across_ci_id, t)
        if ctr == c_ratio_control
            if ctrl_type == pressure_control
                _set_pressure_at_node!(node_across_ci_id, val, ts)
                _set_pressure_at_node!(node_id, cmpr_val * val, ts)
                return
            end
            var1, var2 = _assemble_pipe_contributions_to_node(node_across_ci_id, val, 1/cmpr_val, ts)
            term1 += var1
            term2 += var2
        end
        if ctr == flow_control
            term2 += cmpr_val # inflow positive
        end
    end
    for co in out_c
        ctr, cmpr_val = control(ts, :compressor, co, t)
        node_across_co_id = ref(ts, :compressor, co)["to_node"]
        ctrl_type, val = control(ts, :node, node_across_co_id, t)
        if ctr == c_ratio_control
            if ctrl_type == pressure_control
                _set_pressure_at_node!(node_across_co_id, val, ts)
                _set_pressure_at_node!(node_id, val/cmpr_val, ts)
                return
            end
            var1, var2 = _assemble_pipe_contributions_to_node(node_across_co_id, val, cmpr_val, ts)
            term1 += var1
            term2 += var2
        end
        if ctr == flow_control
            term2 += (-1.0 * cmpr_val) # outflow negative
        end
        if ctr == discharge_pressure_control
            var1, var2 = _assemble_pipe_contributions_to_node(node_across_co_id, val, 1.0, ts)
            discharge_pressure_prev = ref(ts, :node, node_across_co_id, "pressure")
            discharge_rho_prev = get_density(ts, ref(ts, :node, node_across_co_id, "pressure"))
            discharge_rho = get_density(ts, cmpr_val)
            term2 += var1 * (discharge_rho_prev - discharge_rho) + var2
            _set_pressure_at_node!(node_across_co_id, cmpr_val, ts)
        end
    end
    return term1, term2
end

"""
    Assembles all compressor contributions to a node with incoming discharge pressure-controlled compressors
"""
function _assemble_compressor_contributions_to_node_with_incoming_discharge_pressure_control!(compressor_id::Int64, ts::TransientSimulator)::Union{Nothing, Tuple{Real, Real}}
    t = ref(ts, :current_time)
    node_id = ref(ts, :compressor, compressor_id)["fr_node"]
    base_node_id = ref(ts, :compressor, compressor_id)["to_node"]
    ctrl_type, val = control(ts, :node, node_id, t)
    if ctrl_type == pressure_control
        _set_pressure_at_node!(node_id, val, ts)
        # in calling function, check and set pressure at base vertex and across compressors
        return
    end
    out_c = ref(ts, :outgoing_compressors)[base_node_id]
    in_c = ref(ts, :incoming_compressors)[base_node_id]
    term2 = 0.0
    term1 = 0.0
    for ci in in_c
        ctr, cmpr_val = control(ts, :compressor, ci, t)
        if ctr  == c_ratio_control
            node_across_ci_id = ref(ts, :compressor, ci)["fr_node"]
            ctrl_type, val = control(ts, :node, node_across_ci_id, t)
            # if vertex were a slack node, injection is unknown
            @assert ctrl_type == flow_control
            var1, var2 = _assemble_pipe_contributions_to_node(node_across_ci_id, val, 1/cmpr_val, ts)
            term1 += var1
            term2 += var2
        end
        if ctr == flow_control
            term2 += cmpr_val #incoming flow positive
        end
    end
    for co in out_c
        ctr, cmpr_val = control(ts, :compressor, co, t)
        if ctr  == c_ratio_control
            node_across_co_id = ref(ts, :compressor, co)["to_node"]
            ctrl_type, val = control(ts, :node, node_across_co_id, t)
            # if vertex were a slack node, injection is unknown
            @assert ctrl_type == flow_control
            var1, var2 = _assemble_pipe_contributions_to_node(node_across_co_id, val, cmpr_val, ts)
            term1 += var1
            term2 += var2
        end
        if ctr == flow_control
            term2 += (-1.0 * cmpr_val) # outflow
        end
        if ctr == discharge_pressure_control
            var1, var2 = _assemble_pipe_contributions_to_node(node_across_co_id, val, 1.0, ts)
            out_discharge_pressure_prev = ref(ts, :node, node_across_co_id, "pressure")
            out_discharge_rho_prev = ref(ts, :node, node_across_co_id, "density")
            out_discharge_rho = get_density(ts, cmpr_val)
            term2 += var1 * (out_discharge_rho_prev - out_discharge_rho) + var2
        end
    end
    return term1, term2
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
    c = nominal_values(ts, :euler_num) / ( nominal_values(ts, :mach_num) )^2
    beta = ref(ts, :pipe, pipe_id, "friction_factor") / (2 * ref(ts, :pipe, pipe_id, "diameter"))
    a_vec = params(ts, :dt) * beta ./ (rho[2:n] + rho[1:n-1])
    y_vec = phi[2:n] - (c * params(ts, :dt) / ref(ts, :pipe, pipe_id, "dx")) *
            ( get_pressure(ts, rho[2:n]) - get_pressure(ts, rho[1:n-1]) ) - a_vec .* phi[2:n] .* abs.(phi[2:n])
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
    dx = ref(ts, :pipe, pipe_id)["dx"]
    dt = params(ts, :dt)
    n = ref(ts, :pipe, pipe_id)["num_discretization_points"]
    from_node_id = ref(ts, :pipe, pipe_id)["fr_node"]
    rho_prev = get_density(ts, ref(ts, :node, from_node_id)["pressure_previous"])
    rho = get_density(ts, ref(ts, :node, from_node_id)["pressure"])
    # at (n + 1/2) level
    ref(ts, :pipe, pipe_id)["fr_mass_flux"] = ref(ts, :pipe, pipe_id)["fr_minus_mass_flux"] + (rho - rho_prev) * (dx/dt)
    ref(ts, :pipe, pipe_id)["mass_flux_profile"][1] =  ref(ts, :pipe, pipe_id)["fr_mass_flux"]
    # at (n + 1) level
    ref(ts, :pipe, pipe_id)["density_profile"][1] = rho
    to_node_id = ref(ts, :pipe, pipe_id)["to_node"]
    rho_prev = get_density(ts, ref(ts, :node, to_node_id)["pressure_previous"])
    rho = get_density(ts, ref(ts, :node, to_node_id)["pressure"])
    #at (n + 1/2) level
    ref(ts, :pipe, pipe_id)["to_mass_flux"] = ref(ts, :pipe, pipe_id)["to_minus_mass_flux"] + (rho_prev - rho) * (dx/dt)
    ref(ts, :pipe, pipe_id)["mass_flux_profile"][n+1] =  ref(ts, :pipe, pipe_id)["to_mass_flux"]
    # at (n + 1) level
    ref(ts, :pipe, pipe_id)["density_profile"][n] = rho
    return
end

"""
    Reset node flags before time stepping
"""
function _reset_node_flag!(ts::TransientSimulator, node_id::Int64)
    @assert ref(ts, :node, node_id, "is_updated") == true
    ref(ts, :node, node_id)["is_updated"] = false
    return
end
