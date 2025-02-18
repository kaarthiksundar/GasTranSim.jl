function build_ic(data::Dict{String,Any})::Dict{Symbol,Any}
    ic = Dict{Symbol,Any}()

    ic[:node] = Dict() 
    ic[:pipe] = Dict("mass_flow" => Dict(), "pressure" => Dict())
    ic[:compressor] = Dict() 


    
    if haskey(data, "initial_nodal_pressure") && haskey(data, "nodal_pressure")
        (throw(ICException("Both keys nodal pressure and initial_nodal_pressure present")))
    end

    if haskey(data, "initial_nodal_pressure")
        nodal_key = "initial_nodal_pressure"
    else
        nodal_key = "nodal_pressure"
    end

    for (i, value) in get(data, nodal_key, [])
        id = parse(Int64, i)
        ic[:node][id] = value
    end
    


    # this is not a mandatory initial condition - will be computed if not provided
    if haskey(data, "initial_compressor_flow") && haskey(data, "compressor_flow")
        (throw(ICException("Both keys compressor_flow  and initial_compressor_flow present")))
    end

    if haskey(data, "initial_compressor_flow")
        comp_key = "initial_compressor_flow"
    else
        comp_key = "compressor_flow"
    end

    for (i, value) in get(data, comp_key, [])
        id = parse(Int64, i)
        ic[:compressor][id] = value
    end
    


    if haskey(data, "initial_pipe_flow") && haskey(data, "pipe_flow")
        (throw(ICException("Both keys pipe_flow  and initial_pipe_flow present")))
    end

    if haskey(data, "initial_pipe_flow")
        pipe_key = "initial_pipe_flow"
    else
        pipe_key = "pipe_flow"
    end

    for (i, value) in get(data, pipe_key, [])
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

    # compressor flow initial condition is not mandatory
    (isempty(ic[:node])) && (throw(ICException("nodal pressure"))) 
    (isempty(ic[:pipe]["mass_flow"])) && (throw(ICException("pipe flow")))
    
    return ic
end 
