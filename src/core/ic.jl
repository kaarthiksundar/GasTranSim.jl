function build_ic(data::Dict{String,Any})::Dict{Symbol,Any}
    ic = Dict{Symbol,Any}()

    ic[:node] = Dict() 
    ic[:pipe] = Dict("mass_flow" => Dict(), "pressure" => Dict())
    ic[:compressor] = Dict() 

    for (i, value) in get(data, "initial_nodal_pressure", [])
        id = parse(Int64, i)
        ic[:node][id] = value
    end 

    # this is not a mandatory initial condition - will be computed if not provided
    for (i, value) in get(data, "initial_compressor_flow", [])
        id = parse(Int64, i)
        ic[:compressor][id] = value
    end 

    for (i, value) in get(data, "initial_pipe_flow", [])
        id = parse(Int64, i)
        L = data["pipes"][i]["length"]
        if (isa(value, Number))
            distance = [0.0, L]
            val = [value, value]
            ic[:pipe]["mass_flow"][id] = Spline1D(distance, val, k=1)
        else 
            distance = value["distance"]
            val = value["value"]
            ic[:pipe]["mass_flow"][id] = Spline1D(distance, val, k=1)
        end 
    end 

    for (i, value) in get(data, "initial_pipe_pressure", [])
        id = parse(Int64, i)
        distance = value["distance"]
        val = value["value"]
        ic[:pipe]["pressure"][id] = Spline1D(distance, val, k=1)
    end 

    if isempty(get(data, "initial_pipe_pressure", []))
        for (i, pipe) in data["pipes"]
            id = parse(Int64, i)
            L = pipe["length"]
            distance = [0.0, L]
            fr_node = pipe["from_node"]
            to_node = pipe["to_node"]
            fr_pressure = ic[:node][fr_node]
            to_pressure = ic[:node][to_node]
            val = [fr_pressure, to_pressure]
            ic[:pipe]["pressure"][id] = Spline1D(distance, val, k=1)
        end 
    end 

    # compressor flow initial condition is not mandatory
    (isempty(ic[:node])) && (@error "nodal pressure initial condition missing") 
    (isempty(ic[:pipe]["mass_flow"])) && (@error "pipe flow initial condition missing")
    
    return ic
end 
