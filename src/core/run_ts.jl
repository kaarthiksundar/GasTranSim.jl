function run_simulator!(ts::TransientSimulator; 
    run_type::Symbol = :serial, 
    load_adjust::Bool = false,
    showprogress::Bool = true, 
    progress_dt = 1.0)


    minimum_pressure_limit = params(ts, :minimum_pressure_limit)
    ts.params[:load_adjust] = load_adjust
    if  params(ts, :load_adjust) == true && !(minimum_pressure_limit > 0)
        throw(DomainError(minimum_pressure_limit, "load adjustment requires minimum_pressure_limit > 0"))  
    end

    (params(ts, :load_adjust) == true) && (ts.ref[:load_reduction_nodes] = Vector{Int64}())

    output_state = initialize_output_state(ts)
    dt = params(ts, :dt)
    t_f = params(ts, :t_f)
    t_0 = params(ts, :t_0)
    num_steps = Int(round((t_f-t_0)/dt))
    output_data = OutputData(ts)
    prog = Progress(num_steps;
        dt = progress_dt,
        barglyphs = get_barglyphs(),
        barlen = 10, 
        enabled = showprogress, 
        desc = "Sim. progress: ")
    if showprogress == false
        prog = ProgressUnknown(desc="Sim. status ", spinner=true)
    end
    nodal_pressure_previous = form_nodal_pressure_vector(ts)
    for _ in 1:num_steps
    	advance_current_time!(ts, dt)
    	#  if current_time is where some disruption occurs, modify ts.ref now
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
        nodal_pressure_current = form_nodal_pressure_vector(ts)
        error = maximum( abs.(nodal_pressure_current .- nodal_pressure_previous) )/length(nodal_pressure_current)
        
        # commented out temporarily, can be used for debugging code later 
        # if error < 1e-3
        #     println(error)
        # end
    end

    finish!(prog)
    update_output_data!(ts, output_state, output_data)
    populate_solution!(ts, output_data)
end

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

function form_nodal_pressure_vector(ts::TransientSimulator)::Vector{Float64}
    key_array = sort(collect(keys(ref(ts, :node))))
    pressure_vector = Vector{Float64}(undef, length(key_array))
    for i in eachindex(key_array)
        pressure_vector[i] = ref(ts, :node, key_array[i])["pressure"]
    end
    return pressure_vector
end

