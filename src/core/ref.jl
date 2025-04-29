function add_components_to_ref!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

    for (i, node) in get(data, "nodes", [])
        name = :node
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == get(node, "id", get(node, "node_id", []))
        ref[name][id]["id"] = id
        ref[name][id]["is_slack"] = node["slack_bool"]
        ref[name][id]["is_updated"] = false
        ref[name][id]["pressure"] = NaN
        ref[name][id]["pressure_previous"] = NaN
        ref[name][id]["injection"] = NaN
        ref[name][id]["load_reduction"] = 0.0
    end

    for (i, pipe) in get(data, "pipes", [])
        name = :pipe
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == get(pipe, "id", get(pipe, "pipe_id", []))
        ref[name][id]["id"] = id
        ref[name][id]["fr_node"] = get(pipe, "from_node", get(pipe, "fr_node", false))
        if (ref[name][id]["fr_node"] == false)
            throw(MissingDataException("from node for pipe $i "))
        end
        ref[name][id]["to_node"] = pipe["to_node"]
        ref[name][id]["diameter"] = pipe["diameter"]
        ref[name][id]["area"] = pipe["area"]
        ref[name][id]["length"] = pipe["length"]
        ref[name][id]["friction_factor"] = pipe["friction_factor"]
        ref[name][id]["num_discretization_points"] = 0
        ref[name][id]["fr_mass_flux"] = NaN
        ref[name][id]["to_mass_flux"] = NaN
        ref[name][id]["fr_minus_mass_flux"] = NaN
        ref[name][id]["to_minus_mass_flux"] = NaN
    end

    for (i, compressor) in get(data, "compressors", [])
        name = :compressor
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == get(compressor, "id", get(compressor, "comp_id", []))
        ref[name][id]["id"] = id
        ref[name][id]["to_node"] = compressor["to_node"]
        ref[name][id]["fr_node"] =
            get(compressor, "from_node", get(compressor, "fr_node", false))
        if (ref[name][id]["fr_node"] == false)
            throw(MissingDataException("from node for compressor $i "))
        end
        ref[name][id]["control_type"] = unknown_control
        ref[name][id]["c_ratio"] = NaN
        ref[name][id]["discharge_pressure"] = NaN
        ref[name][id]["flow"] = NaN
    end
    return
end

function add_pipe_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

    ref[:incoming_pipes] = Dict{Int64,Vector{Int64}}(i => [] for i in keys(ref[:node]))
    ref[:outgoing_pipes] = Dict{Int64,Vector{Int64}}(i => [] for i in keys(ref[:node]))

    for (id, pipe) in ref[:pipe]
        push!(ref[:incoming_pipes][pipe["to_node"]], id)
        push!(ref[:outgoing_pipes][pipe["fr_node"]], id)
    end
    return
end

function add_compressor_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_compressors] =
        Dict{Int64,Vector{Int64}}(i => [] for i in keys(ref[:node]))
    ref[:outgoing_compressors] =
        Dict{Int64,Vector{Int64}}(i => [] for i in keys(ref[:node]))

    for (id, compressor) in get(ref, :compressor, [])
        push!(ref[:incoming_compressors][compressor["to_node"]], id)
        push!(ref[:outgoing_compressors][compressor["fr_node"]], id)
    end
    return
end

function add_current_time_to_ref!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:current_time] = 0.0
    return
end

function build_ref(data::Dict{String,Any}; ref_extensions = [])::Dict{Symbol,Any}

    ref = Dict{Symbol,Any}()
    add_components_to_ref!(ref, data)

    for extension in ref_extensions
        extension(ref, data)
    end

    return ref
end
