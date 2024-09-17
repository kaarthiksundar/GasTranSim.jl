function run_simulator!(ts::TransientSimulator; 
    run_type::Symbol = :sync, showprogress::Bool = false)
    output_state = initialize_output_state(ts)
    dt = params(ts, :dt)
    t_f = params(ts, :t_f)
    t_0 = params(ts, :t_0)
    num_steps = Int(round((t_f-t_0)/dt))
    output_data = OutputData(ts)
    node_array, _ = order_nodes_for_compressor_flow_computation(ts)
    @showprogress enabled=showprogress desc="Simulation progress: " for _ in 1:num_steps
    	advance_current_time!(ts, dt)
    	#  if current_time is where some disruption occurs, modify ts.ref now
    	advance_pipe_density_internal!(ts, run_type) # (n+1) level
    	advance_node_pressure_mass_flux!(ts, run_type) # pressure (n+1), flux (n+1/2)
    	advance_pipe_mass_flux_internal!(ts, run_type) # (n + 1 + 1/2) level
        calculate_compressor_flows!(ts, node_array)
    	#  if current_time is one where output needs to be saved, check and do now
        update_output_state!(ts, output_state)
    end
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
        (ref(ts, :node, key)["is_level_2"] == true) && (continue)
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

function order_nodes_for_compressor_flow_computation(ts::TransientSimulator)::Tuple{Vector{Int64}, Vector{Int64}}

    node_array = Vector{Int64}()
    for (node_id, _) in ref(ts, :node)
        if abs(ref(ts, :node, node_id)["level"]) == 1 && ref(ts, :node, node_id)["is_slack"] == false
            push!(node_array, node_id)
        end
    end

    all_compressors = collect(keys(ref(ts, :compressor)))
    missing_compressors = Vector{Int64}()

    num_compressors = length(all_compressors)
    cmpr_array = Vector{Int64}()
    node_array_used = Vector{Int64}()

    for node_id in node_array
        
        if ref(ts, :node, node_id)["level"] == 1
            # find the outgoing compressor
            out_c = ref(ts, :outgoing_compressors)[node_id]
            @assert length(out_c) == 1
            if out_c[1] in cmpr_array
                continue
            else
                push!(cmpr_array, out_c[1])
                push!(node_array_used, node_id)
            end
        end
        if ref(ts, :node, node_id)["level"] == -1
            # find the incoming compressor
            in_c = ref(ts, :incoming_compressors)[node_id]
            @assert length(in_c) == 1
            if in_c[1] in cmpr_array
                continue
            else
                push!(cmpr_array, in_c[1])
                push!(node_array_used, node_id)
            end
        end

        if length(cmpr_array) == num_compressors
            break
        end
    end

    if length(cmpr_array) < num_compressors

        @info "The flows in some compressors could not be calculated"
        for id ∈  all_compressors
            if id ∉ cmpr_array
                push!(missing_compressors, id)
            end
        end
    end
    # can use missing_compressors array to calculate using level 2 nodes
    return node_array_used, missing_compressors
end

function calculate_compressor_flows!(ts::TransientSimulator, node_array::Vector)

    for node_id in node_array
        if ref(ts, :node, node_id)["is_slack"] == true
            continue
        end

        t = ref(ts, :current_time)
        ctrl_type, withdrawal = control(ts, :node, node_id, t)
        @assert ctrl_type == flow_control
        _, net_injection = _assemble_pipe_contributions_to_node(node_id, withdrawal, 1.0, ts)

        if ref(ts, :node, node_id)["level"] == 1
            # find the outgoing compressor
            out_c = ref(ts, :outgoing_compressors)[node_id]
            @assert length(out_c) == 1
            ref(ts, :compressor, out_c[1])["flow"] = net_injection
        end
        if ref(ts, :node, node_id)["level"] == -1
            # find the incoming compressor
            in_c = ref(ts, :incoming_compressors)[node_id]
            @assert length(in_c) == 1
            ref(ts, :compressor, in_c[1])["flow"] = -net_injection
        end
    end

    return
end
