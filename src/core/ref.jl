function add_components_to_ref!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

    for (i, node) in get(data, "nodes", [])
        name = :node
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == node["node_id"]
        ref[name][id]["id"] = id
        ref[name][id]["is_slack"] = node["slack_bool"]
        ref[name][id]["is_updated"] = false
        ref[name][id]["pressure"] = NaN 
        ref[name][id]["pressure_previous"] = NaN
        ref[name][id]["injection"] = NaN
    end

    for (i, pipe) in get(data, "pipes", [])
        name = :pipe
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == pipe["pipe_id"]
        ref[name][id]["id"] = id
        ref[name][id]["fr_node"] = pipe["from_node"]
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
        @assert id == compressor["comp_id"]
        ref[name][id]["id"] = id
        ref[name][id]["to_node"] = compressor["to_node"]
        ref[name][id]["fr_node"] = compressor["from_node"]
        ref[name][id]["control_type"] = unknown_control
        ref[name][id]["c_ratio"] = NaN
        ref[name][id]["discharge_pressure"] = NaN
        ref[name][id]["flow"] = NaN
    end
    return
end

function add_pipe_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

    ref[:incoming_pipes] = Dict{Int64, Vector{Int64}}()
    ref[:outgoing_pipes] = Dict{Int64, Vector{Int64}}()

    for (i, val) in ref[:node]
        ref[:incoming_pipes][i] = []
        ref[:outgoing_pipes][i] = []
    end

    for (id, pipe) in ref[:pipe]
        push!(ref[:incoming_pipes][pipe["to_node"]], id)
        push!(ref[:outgoing_pipes][pipe["fr_node"]], id)
    end
    return
end

function add_compressor_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

    ref[:incoming_compressors] = Dict{Int64, Vector{Int64}}()
    ref[:outgoing_compressors] = Dict{Int64, Vector{Int64}}()
    
    for (i, val) in ref[:node]
        ref[:incoming_compressors][i] = []
        ref[:outgoing_compressors][i] = []
    end

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

function build_ref(data::Dict{String,Any};
    ref_extensions=[])::Dict{Symbol,Any}

    ref = Dict{Symbol,Any}()
    add_components_to_ref!(ref, data)

    for extension in ref_extensions
        extension(ref, data)
    end

    return ref
end
