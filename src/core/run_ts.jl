function run_simulator!(ts::TransientSimulator; run_type = :sync)
    out_int = initialize_output_struc(ts)
    dt = params(ts, :dt)
    while ts.ref[:current_time] < 4.0 * dt
    	advance_current_time!(ts, dt)
    	#  if current_time is where some disruption occurs, modify ts.ref now
    	advance_pipe_density_internal!(ts, run_type) #(n+1) level
    	advance_node_pressure_mass_flux!(ts, run_type) #pressure (n+1), flux (n+1/2)
    	advance_mass_flux_internal!(ts, run_type) # (n+1 + 1/2) level
    	#  if current_time is one where output needs to be saved, check and do now
    	update_output_struc!(ts, out_int)
    	 err_arr = []
    	 for key in collect(keys(ts.ref[:node]))
    	 	push!(err_arr, abs(ts.ref[:node][key]["pressure_previous"] - ts.ref[:node][key]["pressure"]) )
    	 end
    	 @show maximum(err_arr)
    end
    # flux profile, density profiles saved to restart time marching
    #out = create_output(ts, out_int)
end

function advance_current_time!(ts::TransientSimulator, tau::Type{<:Real})
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
    for (key, junction) in ref(ts, :node)
        (ref(ts, :node, key)["is_updated"] == true) && (continue)
        # p(t), but q(t - dt/2) taken care of inside
        ctrl_type, val = control(ts, :node, key, t)
        if ctrl_type == pressure_control
            _set_pressure_at_node!(key, val, ts)
            _set_pressure_at_node_across_compressors!(key, val, ts)
        elseif ctrl_type == flow_control
            _solve_for_pressure_at_node_and_neighbours!(key, val, ts)
        else
            @error "control type unknown at advance_pressure_mass_flux_node!"
        end
    end

    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_compute_end_fluxes_densities!, ts, key_array, run_type)
    key_array = collect(keys(ref(ts, :node)))
    _execute_task!(_reset_vertex_flag!, ts, key_array, run_type)
    return
end

function advance_mass_flux_internal!(ts::TransientSimulator, run_type::Symbol)
    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_advance_pipe_mass_flux_internal!, ts, key_array, run_type)
    return
end
