function run_simulator!(ts::TransientSimulator; 
    run_type::Symbol = :sync, 
    showprogress::Bool = false, 
    progress_dt = 1.0)
    output_state = initialize_output_state(ts)
    dt = params(ts, :dt)
    t_f = params(ts, :t_f)
    t_0 = params(ts, :t_0)
    num_steps = Int(round((t_f-t_0)/dt))
    output_data = OutputData(ts)
    prog = Progress(num_steps;
        dt = progress_dt,
        barlen = 10, 
        enabled = showprogress, 
        desc = "Sim. progress: ")
    if showprogress == false
        prog = ProgressUnknown(desc="Sim. status", spinner=true)
    end 
    for _ in 1:num_steps
    	advance_current_time!(ts, dt)
    	#  if current_time is where some disruption occurs, modify ts.ref now
        advance_nodal_pressure!(ts)

        # advance_nodal_pressure!() #(n+1) level
          # go to each slack node and update pressure, mark update node flag
            # advance_compressor_nodes!()  #(n+1) level
            # if  slack node has compressor, update other end vertex, mark update node flag
            # if any compressor has discharge pressure set, use to update vertex, mark update node flag
          # now go to other nodes
            #calculate density based on conservation eq with chosen fictitious vertex vol and edge mass flux    from previous step # check - does this need theta or smhtg to be second order accurate in time  and space to match pipe scheme ?
          # 
        # 
        # advance_pipes!() #(n+1) level
            # global id for pressure, mass flux
            # assemble sparse mat, load vec
            # solve
            # postprocess and fill in nodal pressure, density, fr_flux, to_flux
        # fill in compressor fluuxes, calculate compressor vals by 
            # use vol of vertex to now find compressor flows
        
        #=
    	advance_pipe_density_internal!(ts, run_type) # (n+1) level
    	advance_node_pressure_mass_flux!(ts, run_type) # pressure (n+1), flux (n+1/2)
    	advance_pipe_mass_flux_internal!(ts, run_type) # (n + 1 + 1/2) level
        _compute_compressor_flows!(ts)
    	#  if current_time is one where output needs to be saved, check and do now
        update_output_state!(ts, output_state)
        if showprogress == false
            next!(prog, spinner="🌑🌒🌓🌔🌕🌖🌗🌘")
        else 
            next!(prog)
        end 
    end
    finish!(prog)
    update_output_data!(ts, output_state, output_data)
    populate_solution!(ts, output_data)
        =#
    end
    # update_output_data!(ts, output_state, output_data)
    # populate_solution!(ts, output_data)
end


function advance_nodal_pressure!(ts::TransientSimulator)
    # DO NOT parallelize this (race condition)
    for (key, junction) in ref(ts, :node)
        (ref(ts, :node, key)["is_updated"] == true) && (continue)
        # (ref(ts, :node, key)["is_level_2"] == true) && (continue)
        # p(t), but q(t) taken care of inside
        ctrl_type, val = control(ts, :node, key, t)
        if ctrl_type == pressure_control
            _set_pressure_at_node!(key, val, ts)
            _set_pressure_at_node_across_compressors!(key, val, ts)
        elseif ctrl_type == flow_control
            _solve_for_pressure_at_node_and_neighbours!(key, val, ts)
        else
            @error "control type unknown at advance_nodal_pressure!"
        end
    end

    return
end
#---------------
function advance_current_time!(ts::TransientSimulator, tau::Real)
    ts.ref[:current_time] += tau
    return
end

function advance_pipe_density_internal!(ts::TransientSimulator, run_type::Symbol)
    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_advance_pipe_density_internal!, ts, key_array, run_type)
    return
end

function advance_node_pressure_mass_flux!(ts::TransientSimulator, run_type::Symbol)
    t = ref(ts, :current_time)
    # DO NOT parallelize this (race condition)
    # for (key, junction) in ref(ts, :node)
    for (key, junction) in ref(ts, :node)
        (ref(ts, :node, key)["is_updated"] == true) && (continue)
        (ref(ts, :node, key)["level"] == 2) && (continue)
        # p(t), but q(t - dt/2) taken care of inside
        ctrl_type, val = control(ts, :node, key, t)
        if ctrl_type == pressure_control
            _set_pressure_at_node!(key, val, ts)
            _set_pressure_at_node_across_compressors!(key, val, ts)
        elseif ctrl_type == flow_control
            _solve_for_pressure_at_node_and_neighbours!(key, val, ts)
        else
            throw(ControlException("control type unknown at advance_pressure_mass_flux_node"))
        end
    end

    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_compute_pipe_end_fluxes_densities!, ts, key_array, run_type)
    key_array = collect(keys(ref(ts, :node)))
    _execute_task!(_reset_node_flag!, ts, key_array, run_type)
    return
end

function advance_pipe_mass_flux_internal!(ts::TransientSimulator, run_type::Symbol)
    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_advance_pipe_mass_flux_internal!, ts, key_array, run_type)
    return
end

