
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

"""
    Sets the pressure and density at a given node to a value
"""
function _set_pressure_at_node!(node_id, pressure, ts::TransientSimulator)
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
function _set_pressure_at_node_across_compressors!(node_id, pressure, ts::TransientSimulator)
    t = ref(ts, :current_time)
    for ci in ref(ts,:incoming_compressors, node_id)
        ctrl_type, val = control(ts, :compressor, ci, t)
        i = ref(ts, :compressor, ci, "fr_node")
        if ctrl_type == c_ratio_control
            _set_pressure_at_node!(i, pressure/val, ts)
        end
    end
    for ci in ref(ts,:outgoing_compressors, node_id)
        ctrl_type, val = control(ts, :compressor, ci, t)
        i = ref(ts, :compressor, ci, "to_node")
        if ctrl_type == c_ratio_control
            _set_pressure_at_node!(i, pressure * val, ts)
        elseif ctrl_type == discharge_pressure_control
            _set_pressure_at_node!(i, val, ts)
        end
    end
    return
end

"""
    Solves for the pressure at the given node (with flow control) and its neighbors
"""
function _solve_for_pressure_at_node_and_neighbours!(node_id, withdrawal, ts::TransientSimulator)
    t = ref(ts, :current_time)
    in_c = ref(ts, :incoming_compressors, node_id)
    index = findfirst(ci -> control(ts, :compressor, ci, t)[1] == discharge_pressure_control, in_c)
    # we know there can be at most one such compressor
    if isa(index, Nothing)
        _calculate_pressure_for_node_without_incoming_discharge_pressure_control!(node_id, withdrawal, ts)
    else
        _calculate_pressure_for_node_with_incoming_discharge_pressure_control!(in_c[index], ts)
    end
    return
end


function _invert_quadratic(a::Type{<:Real}, y::Type{<:Real})::Type{<:Real}
    # can also write as 2 * y / (1 + sqrt( 1 + 4 * a * abs(y) ) )
    return sign(y) * ( -1.0 + sqrt( 1.0 + 4.0 * a * abs(y) ) ) / (2.0 * a)
end

function _advance_pipe_mass_flux_internal!(ts::TransientSimulator, key::Int64)
        rho = ts.ref[:pipe][key]["density_profile"]
        phi = ts.ref[:pipe][key]["mass_flux_profile"]
        n = ts.ref[:pipe][key]["num_discretization_points"]
        beta = ts.ref[:pipe][key]["friction_factor"] / (2 * ts.ref[:pipe][key]["diameter"])

        # for i = 2:n
        #     # even if denom of a is effectively 0, this code does not give error
        #     a = ts.params[:dt] * beta / (rho[i] + rho[i-1])
        #     y = phi[i] - ( ts.params[:dt] / ts.ref[:pipe][key]["dx"] ) *
        #     (get_pressure(ts, rho[i]) - get_pressure(ts, rho[i-1])) - a * phi[i] * abs(phi[i])
        #     phi[i] = _invert_quadratic(a, y)
        # end
        # can vectorize this as

        a_vec = ts.params[:dt] * beta  ./ (rho[2:n] + rho[1:n-1])
        y_vec = phi[2:n] - ( ts.params[:dt] / ts.ref[:pipe][key]["dx"] ) *
            (get_pressure(ts, rho[2:n]) - get_pressure(ts, rho[1:n-1])) - a_vec .* phi[2:n] .* abs.(phi[2:n])
        phi[2:n] = _invert_quadratic.(a_vec, y_vec)

        # update field
        ts.ref[:pipe][key]["fr_minus_mass_flux"] = phi[2]
        ts.ref[:pipe][key]["to_minus_mass_flux"] = phi[n]
        return
end



function _calculate_pressure_for_node_without_incoming_discharge_pressure_control!(vertex_key::Int64, val::Float64, ts::TransientSimulator)

    s  = _assemble_vertex_contribution_compressors_without_incoming_discharge_pressure_control!(vertex_key, ts)
    if isa(s, Nothing)
        # means vertex pressure is already set
        prs_val = ts.ref[:node][vertex_key]["pressure"]
        _set_pressure_at_node_across_compressors!(vertex_key, prs_val, ts)
        return
    end
    s1, s2 = s
    t1, t2 = _assemble_vertex_contributions_pipes(vertex_key, val, 1.0, ts)


    rho_prev = get_density(ts, ts.ref[:node][vertex_key]["pressure"])
    rho = rho_prev + (t2 + s2)/(t1 + s1)
    prs_val = get_pressure(ts, rho)

    _set_pressure_at_node!(vertex_key, prs_val, ts)
    _set_pressure_at_node_across_compressors!(vertex_key, prs_val, ts)

    return

end




function _assemble_vertex_contributions_pipes(vertex_key::Int64, val::Float64, mult_factor::Float64, ts::TransientSimulator)::Tuple{Real, Real}

    out_p = ts.ref[:outgoing_pipes][vertex_key]
    in_p = ts.ref[:incoming_pipes][vertex_key]
    term1 = 0.0
    term2 = -1.0 * val # in input data, withdrawal is positive, but we want inflow positive
    pipe = ts.ref[:pipe]

    for p in out_p
        term1 += mult_factor * pipe[p]["dx"] * pipe[p]["area"] / ts.params[:dt]
        term2 -= pipe[p]["fr_minus_mass_flux"] * pipe[p]["area"]
    end

    for p in in_p
        term1 += mult_factor * pipe[p]["dx"] * pipe[p]["area"] / ts.params[:dt]
        term2 += pipe[p]["to_minus_mass_flux"] * pipe[p]["area"]
    end



    return term1, term2
end


function _assemble_vertex_contribution_compressors_without_incoming_discharge_pressure_control!(vertex::Int64, ts::TransientSimulator)::Union{Nothing, Tuple{Real, Real}}

    out_c = ts.ref[:outgoing_compressors][vertex]
    in_c = ts.ref[:incoming_compressors][vertex]
    term1 = 0.0
    term2 = 0.0
    p_prev = ts.ref[:node][vertex]["pressure"]
    rho_prev = get_density(ts, p_prev)



    for ci in in_c
            ctr, cmpr_val = control(ts, :compressor, ci, ts.ref[:current_time])
            vertex_across_ci = ts.ref[:compressor][ci]["fr_node"]

            if ctr  == c_ratio_control

                ctrl_type, val = control(ts, :node, vertex_across_ci, ts.ref[:current_time])

                if ctrl_type == pressure_control
                    _set_pressure_at_node!(vertex_across_ci, val, ts)
                    _set_pressure_at_node!(vertex, cmpr_val*val, ts)
                    return
                end


                var1, var2 = _assemble_vertex_contribution_pipes(vertex_across_ci, val, 1/cmpr_val, ts)
                term1 += var1
                term2 += var2

            end

            if ctr == flow_control
                term2 += cmpr_val #inflow positive
            end

    end


    for co in out_c
            ctr, cmpr_val = control(ts, :compressor, co, ts.ref[:current_time])
            vertex_across_co = ts.ref[:compressor][co]["to_node"]
            ctrl_type, val = control(ts, :node, vertex_across_co, ts.ref[:current_time])

            if ctr  == c_ratio_control

                if ctrl_type == pressure_control
                    _set_pressure_at_node!(vertex_across_co, val, ts)
                    _set_pressure_at_node!(vertex, val/cmpr_val, ts)
                    return
                end

                var1, var2 = _assemble_vertex_contributions_pipes(vertex_across_co, val, cmpr_val, ts)
                term1 += var1
                term2 += var2

            end

            if ctr == flow_control
                term2 += -1.0 * cmpr_val #outflow negative
            end

            if ctr == discharge_pressure_control
                var1, var2 = _assemble_vertex_contributions_pipes(vertex_across_co, val, 1.0, ts)
                pressure_prev_del = ts.ref[:node][vertex_across_co]["pressure"]
                rho_prev_del = get_density(ts, pressure_prev_del)
                rho_del = get_density(ts, cmpr_val)
                term2 += var1 *(rho_prev_del - rho_del) + var2
                _set_pressure_at_node!(vertex_across_co, cmpr_val, ts)

            end



    end

    return term1, term2

end



function _assemble_vertex_contribution_compressors_with_incoming_discharge_pressure_control!(compressor_id::Int64, ts::TransientSimulator)::Union{Nothing, Tuple{Real, Real}}

    vertex = ts.ref[:compressor][compressor_id]["fr_node"]

    base_vertex = ts.ref[:compressor][compressor_id]["to_node"]

    ctrl_type, val = control(ts, :node, vertex, ts.ref[:current_time])

    if ctrl_type == pressure_control
        _set_pressure_at_node!(vertex, val, ts)
        #in calling function, check and set pressure at base vertex and across compressors
        return
    end

    out_c = ts.ref[:outgoing_compressors][base_vertex]
    in_c = ts.ref[:incoming_compressors][base_vertex]
    term2 = 0.0
    term1 = 0.0

    for ci in in_c
            ctr, cmpr_val = control(ts, :compressor, ci, ts.ref[:current_time])

            if ctr  == c_ratio_control
                vertex_across_ci = ts.ref[:compressor][ci]["fr_node"]
                ctrl_type, val = control(ts, :node, vertex_across_ci, ts.ref[:current_time])

                # if vertex were a slack node, injection is unknown
                @assert ctrl_type == flow_control

                var1, var2 = _assemble_vertex_contribution_pipes(vertex_across_ci, val, 1/cmpr_val, ts)
                term1 += var1
                term2 += var2
            end

            if ctr == flow_control
                term2 += cmpr_val #incoming flow positive
            end

    end


    for co in out_c
            ctr, cmpr_val = control(ts, :compressor, co, ts.ref[:current_time])
            if ctr  == c_ratio_control
                vertex_across_co = ts.ref[:compressor][co]["to_node"]
                ctrl_type, val = control(ts, :node, vertex_across_co, ts.ref[:current_time])

                # if vertex were a slack node, injection is unknown
                @assert ctrl_type == flow_control

                var1, var2 = _assemble_coefficients_pipes(vertex_across_co, val, cmpr_val, ts)
                term1 += var1
                term2 += var2
            end

            if ctr == flow_control
                term2 += -1.0 * cmpr_val #outflow
            end

            if ctr == discharge_pressure_control
                var1, var2 = _assemble_coefficients_pipes(vertex_across_co, val, 1.0, ts)
                out_pressure_prev_del = ts.ref[:node][vertex_across_co]["pressure"]
                out_rho_prev_del = get_density(ts, out_pressure_prev_del)
                out_rho_del = get_density(ts, cmpr_val)
                term2 += var1 *(out_rho_prev_del - out_rho_del) + var2

            end



    end

    return term1, term2

end


function _calculate_pressure_for_node_with_incoming_discharge_pressure_control!(compressor_id::Int64, ts::TransientSimulator)

    vertex = ts.ref[:compressor][compressor_id]["fr_node"]
    base_vertex = ts.ref[:compressor][compressor_id]["to_node"]
    ctrl, pressure_del = control(ts, :compressor, compressor_id, ts.ref[:current_time])
    @assert ctrl == discharge_pressure_control

    s = _assemble_vertex_contribution_compressors_with_incoming_discharge_pressure_control!(compressor_id, ts)
    if isa(s, Nothing)
        #pressure at from end of discharge compressor already set
        _set_pressure_at_node!(base_vertex, pressure_del, ts)
        _set_pressure_at_node_across_compressors!(base_vertex, pressure_del, ts)
        return
    end
    s1, s2 = s

    ctrl_type, val1 = control(ts, :node, vertex, ts.ref[:current_time])
    @assert ctrl_type == flow_control
    t1, t2 = _assemble_vertex_contributions_pipes(vertex, val1, 1.0, ts)

    ctrl_type, val2 = control(ts, :node, base_vertex, ts.ref[:current_time])
    @assert ctrl_type == flow_control
    r1, r2 = _assemble_vertex_contributions_pipes(base_vertex, val2, 1.0, ts)


    rho_del = get_density(ts, pressure_del)
    pressure_prev_del = ts.ref[:node][base_vertex]["pressure"]
    rho_prev_del = get_density(ts, pressure_prev_del)

    rho_prev = get_density(ts, ts.ref[:node][vertex]["pressure"])
    rho = rho_prev + (t2 +  (r1 + s1) * (rho_prev_del - rho_del) + r2 + s2) / t1

    _set_pressure_at_node!(vertex, get_pressure(ts, rho), ts)
    _set_pressure_at_node!(base_vertex, pressure_del, ts)
    _set_pressure_at_node_across_compressors!(base_vertex, pressure_del, ts)
    return

end

function _compute_end_fluxes_densities!(ts::TransientSimulator, pipe_id::Int64)

    dx = ts.ref[:pipe][pipe_id]["dx"]
    dt = ts.params[:dt]
    n = ts.ref[:pipe][pipe_id]["num_discretization_points"]

    from_vertex = ts.ref[:pipe][pipe_id]["fr_node"]
    rho_prev = get_density(ts, ts.ref[:node][from_vertex]["pressure_previous"])
    rho = get_density(ts, ts.ref[:node][from_vertex]["pressure"])
    # at n+ 1/2 level
    ts.ref[:pipe][pipe_id]["fr_mass_flux"] =  ts.ref[:pipe][pipe_id]["fr_minus_mass_flux"] + (rho - rho_prev) * (dx/dt)
    # at n+1 level
    ts.ref[:pipe][pipe_id]["density_profile"][1] = rho

    to_vertex = ts.ref[:pipe][pipe_id]["to_node"]
    rho_prev = get_density(ts, ts.ref[:node][to_vertex]["pressure_previous"])
    rho = get_density(ts, ts.ref[:node][to_vertex]["pressure"])
    #at n+ 1/2 level
    ts.ref[:pipe][pipe_id]["to_mass_flux"] = ts.ref[:pipe][pipe_id]["to_mass_flux"] + (rho_prev - rho) * (dx/dt)
    # at n+1 level
    ts.ref[:pipe][pipe_id]["density_profile"][n] = rho
    return

end

function _reset_vertex_flag!(ts::TransientSimulator, key::Int64)
    @assert ts.ref[:node][key]["is_updated"] = true # ensures all vertices traversed
    ts.ref[:node][key]["is_updated"] = false
    return
end
