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
        ref[name][id]["area"] = pi * pipe["diameter"] * pipe["diameter"] * 0.25
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
        ref[name][id]["fr_node"] = compressor["from_node"]
        ref[name][id]["to_node"] = compressor["to_node"]
        ref[name][id]["control_type"] = unknown
        ref[name][id]["c_ratio"] = NaN 
        ref[name][id]["discharge_pressure"] = NaN 
        ref[name][id]["flow"] = NaN 
    end 
end 

function add_initial_conditions_to_ref!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

    for (i, value) in get(data, "initial_nodal_pressure", [])
        id = parse(Int64, i)
        @assert id in keys(ref[:node])
        ref[:node][id]["initial_pressure"] = value
    end 

    for (i, value) in get(data, "initial_nodal_flow", [])
        id = parse(Int64, i)
        @assert id in keys(ref[:node])
        ref[:node][id]["initial_injection"] = (-1.0) * value
    end 

    for (i, value) in get(data, "initial_pipe_flow", [])
        id = parse(Int64, i)
        @assert id in keys(ref[:pipe])
        ref[:pipe][id]["initial_mass_flux"] = value / ref[:pipe][id]["area"]
    end 

    for (i, value) in get(data, "initial_pipe_pressure_in", [])
        id = parse(Int64, i)
        @assert id in keys(ref[:pipe])
        ref[:pipe][id]["initial_fr_pressure"] = value 
        
    end 

    for (i, value) in get(data, "initial_pipe_pressure_out", [])
        id = parse(Int64, i)
        @assert id in keys(ref[:pipe])
        ref[:pipe][id]["initial_to_pressure"] = value 
    end 

    for (i, value) in get(data, "initial_compressor_flow", [])
        id = parse(Int64, i)
        @assert id in keys(ref[:compressor])
        ref[:compressor][id]["initial_mass_flow"] = value 
    end 

    for (i, value) in get(data, "initial_compressor_pressure_in", [])
        id = parse(Int64, i)
        @assert id in keys(ref[:compressor])
        ref[:compressor][id]["initial_fr_pressure"] = value 
        
    end 

    for (i, value) in get(data, "initial_compressor_pressure_out", [])
        id = parse(Int64, i)
        @assert id in keys(ref[:compressor])
        ref[:compressor][id]["initial_to_pressure"] = value 
    end
    
    for (i, value) in get(data, "initial_compressor_ratio", [])
        id = parse(Int64, i)
        @assert id in keys(ref[:compressor])
        ref[:compressor][id]["initial_c_ratio"] = value 
    end
end 

function build_ref(data::Dict{String,Any};
    ref_extensions=[])::Dict{Symbol,Any}

    ref = Dict{Symbol,Any}()
    add_components_to_ref!(ref, data) 
    add_initial_conditions_to_ref!(ref, data)
    
    for extension in ref_extensions 
        extension(ref, data)
    end 

    return ref
end 