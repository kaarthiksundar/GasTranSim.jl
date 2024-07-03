function get_data(data_folder::AbstractString; 
    case_name::AbstractString="", 
    case_types::Vector{Symbol}=Symbol[])::Dict{String,Any}
    return parse_data(data_folder; case_name=case_name, case_types=case_types)
end

function initialize_simulator(data_folder::AbstractString;
    case_name::AbstractString="", 
    case_types::Vector{Symbol}=Symbol[],
    kwargs...)::TransientSimulator
    data = parse_data(data_folder; case_name=case_name, case_types=case_types)
    return initialize_simulator(data; kwargs...)
end

function initialize_simulator(data::Dict{String,Any}; eos::Symbol=:ideal)::TransientSimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    ref = build_ref(data, ref_extensions= [
        add_pipe_info_at_nodes!,
        add_compressor_info_at_nodes!,
        add_current_time_to_ref!]
    )
    ref[:current_time] = params[:t_0]

    ts = TransientSimulator(data,
        ref,
        initialize_solution(data),
        nominal_values,
        params,
        build_ic(data),
        build_bc(data),
        get_eos(eos)...
    )

    add_pipe_grid_to_ref!(ts)
    add_node_level_flag!(ts)
    add_node_ordering_for_compressor_flow_computation!(ts)
    initialize_nodal_state!(ts)
    initialize_pipe_state!(ts)
    initialize_compressor_state!(ts)
    return ts
end


function add_pipe_grid_to_ref!(ts::TransientSimulator)
    for (key, pipe) in ref(ts, :pipe)

        dx  = ts.params[:dx]
        ref(ts, :pipe, key)["num_discretization_points"] =   Int(floor(1 + pipe["length"] / dx))
        n = ref(ts, :pipe, key)["num_discretization_points"]
        ref(ts, :pipe, key)["density_profile"] = zeros(Float64, n)
        ref(ts, :pipe, key)["mass_flux_profile"] = zeros(Float64, n)
    end
    return
end

function _evaluate_level_of_node!(ts::TransientSimulator, node_id::Int64)
    if !haskey(ref(ts), :compressor)
        ref(ts, :node, node_id)["level"] = 0 
        return
    end
    if length(ref(ts, :incoming_compressors, node_id)) + length(ref(ts, :outgoing_compressors, node_id)) == 0
        ref(ts, :node, node_id)["level"] = 0 
        return
    end
    if length(ref(ts, :incoming_compressors, node_id)) + length(ref(ts, :outgoing_compressors, node_id)) == 1
        if length(ref(ts, :incoming_compressors, node_id)) == 1
            ref(ts, :node, node_id)["level"] = -1
        end

        if length(ref(ts, :outgoing_compressors, node_id)) == 1
            ref(ts, :node, node_id)["level"] =  1
        end

        return
    end
    for ci in ref(ts, :incoming_compressors, node_id)
        node_across_ci = ref(ts, :compressor, ci, "fr_node")
        if length(ref(ts, :incoming_compressors, node_across_ci)) + length(ref(ts, :outgoing_compressors, node_across_ci))  >  1
            throw(NetworkException("3 compressors are in series"))
        end
    end
    for ci in ref(ts, :outgoing_compressors, node_id)
        node_across_ci = ref(ts, :compressor, ci, "to_node")
        if length(ref(ts, :incoming_compressors, node_across_ci)) + length(ref(ts, :outgoing_compressors, node_across_ci))  > 1
            throw(NetworkException("3 compressors are in series"))
        end
    end
    ref(ts, :node, node_id)["level"] = 2
    return
end

function add_node_level_flag!(ts::TransientSimulator)
    for (node_id, _) in ref(ts, :node)
        _evaluate_level_of_node!(ts, node_id)
    end
    return
end

function add_node_ordering_for_compressor_flow_computation!(ts::TransientSimulator)
    candidate_node_ids = Vector{Int64}()
    for (node_id, _) in ref(ts, :node)
        if abs(ref(ts, :node, node_id)["level"]) == 1 && ref(ts, :node, node_id)["is_slack"] == false
            push!(candidate_node_ids, node_id)
        end
    end

    compressors = get(ref(ts), :compressor, Dict())
    (isempty(compressors)) && (return)
    compressor_ids = keys(compressors) |> collect
    compressor_ids_with_uncomputable_flow = Vector{Int64}()
    num_compressors = length(compressors)
    compressors_accounted_for = Vector{Int64}()
    used_node_ids = Vector{Int64}()

    for node_id in candidate_node_ids
        if ref(ts, :node, node_id)["level"] == 1
            # find the outgoing compressor
            out_c = ref(ts, :outgoing_compressors)[node_id]
            # @assert length(out_c) == 1
            (out_c[1] in compressors_accounted_for) && (continue)
            push!(compressors_accounted_for, out_c[1])
            push!(used_node_ids, node_id)
        end
        if ref(ts, :node, node_id)["level"] == -1
            # find the incoming compressor
            in_c = ref(ts, :incoming_compressors)[node_id]
            # @assert length(in_c) == 1
            (in_c[1] in compressors_accounted_for) && (continue)
            push!(compressors_accounted_for, in_c[1])
            push!(used_node_ids, node_id)
        end
    end 

    # can use missing_compressors array to calculate using level 2 nodes (TODO)
    if length(compressors_accounted_for) < num_compressors
        @info "The flows in some compressors cannot be calculated"
        for id in compressor_ids 
            !(id in compressors_accounted_for) && (push!(compressor_ids_with_uncomputable_flow, id))
        end 
    end 

    ref(ts)[:node_ordering_for_compressor_flow_calculation] = used_node_ids 
    ref(ts)[:compressor_ids_with_uncomputable_flow] = compressor_ids_with_uncomputable_flow
end 

function initialize_nodal_state!(ts::TransientSimulator)
    for (key, _) in ref(ts, :node)
        pressure = initial_nodal_pressure(ts, key)
        ref(ts, :node, key)["pressure"] = pressure
    end
    return
end


function initialize_pipe_state!(ts::TransientSimulator)
    is_steady = false 
    if isempty(ts.initial_conditions[:pipe]["pressure"])
        is_steady = true
    end 
    for (key, pipe) in ref(ts, :pipe)
        area = pipe["area"]
        n = pipe["num_discretization_points"]
        dx = ts.params[:dx]
        L = pipe["length"]
        fr_node = pipe["fr_node"]
        to_node = pipe["to_node"]
        if is_steady
            initial_mass_flux = initial_pipe_mass_flow(ts, key)(0.0) / area
            fill!(pipe["mass_flux_profile"], initial_mass_flux)
            initial_fr_pressure = ref(ts, :node, fr_node, "pressure")
            initial_to_pressure = ref(ts, :node, to_node, "pressure")
            density_at_first_sq = get_density(ts, initial_fr_pressure) ^ 2
            density_at_last_sq = get_density(ts, initial_to_pressure) ^ 2
            dL = dx / L
            pipe["density_profile"][1:n] = [
                sqrt(
                    density_at_last_sq * (i  - 1) * dL +
                    density_at_first_sq * (n  - i) * dL
                ) for i in 1:n
            ]
            # pipe["fr_minus_mass_flux"] = initial_mass_flux # (dx/2)
            # pipe["to_minus_mass_flux"] = initial_mass_flux # L-(dx/2)
            pipe["fr_mass_flux"] = initial_mass_flux # -(dx/2)
            pipe["to_mass_flux"] = initial_mass_flux # L+(dx/2)
        else 
            flow_spl = initial_pipe_mass_flow(ts, key)
            pressure_spl = initial_pipe_pressure(ts, key)
            x_rho = LinRange(0, L, n)
            # x_mid = x_rho[1:n-1] .+ dx/2.0  # modify since no longer staggered
            get_coeffs(flow_spl)[1]
            pipe["mass_flux_profile"] = [
                get_coeffs(flow_spl)[1], 
                [flow_spl(x) for x in x_rho]..., 
                get_coeffs(flow_spl)[end]
            ] ./ area
            # pipe["fr_minus_mass_flux"] = pipe["mass_flux_profile"][2]
            # pipe["to_minus_mass_flux"] = pipe["mass_flux_profile"][end-1]
            pipe["fr_mass_flux"] = pipe["mass_flux_profile"][1]
            pipe["to_mass_flux"] = pipe["mass_flux_profile"][end]
            pipe["density_profile"][1:n] = [
                get_density(ts, pressure_spl(x)) for x in x_rho
            ]
        end 
    end
    return
end

function _compute_compressor_flows!(ts::TransientSimulator)
    for node_id in get(ref(ts), :node_ordering_for_compressor_flow_calculation, [])
        t = ref(ts, :current_time)
        ctrl_type, withdrawal = control(ts, :node, node_id, t)
        @assert ctrl_type == flow_control
        _, net_injection = _assemble_pipe_contributions_to_node(node_id, withdrawal, 1.0, ts)

        if ref(ts, :node, node_id)["level"] == 1
            # find the outgoing compressor
            out_c = ref(ts, :outgoing_compressors)[node_id]
            # @assert length(out_c) == 1
            ref(ts, :compressor, out_c[1])["flow"] = net_injection
        end
        if ref(ts, :node, node_id)["level"] == -1
            # find the incoming compressor
            in_c = ref(ts, :incoming_compressors)[node_id]
            # @assert length(in_c) == 1
            ref(ts, :compressor, in_c[1])["flow"] = -net_injection
        end
    end
    return
end

function initialize_compressor_state!(ts::TransientSimulator)
    isempty(get(ref(ts), :compressor, [])) && return
    if has_compressor_initial_flows(ts)
        for (key, _) in ref(ts, :compressor)
            flow = initial_compressor_flow(ts, key)
            ref(ts, :compressor, key)["flow"] = flow
        end
        return
    end  
    @info "compressor does not have initial flows, computing them"
    _compute_compressor_flows!(ts)
    return
end