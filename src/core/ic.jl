function build_ic(data::Dict{String,Any})::Dict{Symbol,Any}
    ic = Dict{Symbol,Any}()

    ic[:node] = Dict() 
    ic[:pipe] = Dict("mass_flow" => Dict(), "pressure" => Dict())
    ic[:compressor] = Dict() 

    # initial conditions file should not have both "initial_nodal_pressure" and "nodal_pressure"
    if haskey(data, "initial_nodal_pressure") && haskey(data, "nodal_pressure")
        (throw(ICException("both keys nodal pressure and initial_nodal_pressure present")))
    end

    nodal_pressure_key = haskey(data, "initial_nodal_pressure") ? "initial_nodal_pressure" : "nodal_pressure"
    
    for (i, value) in get(data, nodal_pressure_key, [])
        id = parse(Int64, i)
        ic[:node][id] = value
    end
    
    # this is not a mandatory initial condition - will be computed if not provided
    # initial conditions file should not have both "compressor_flow" and "initial_compressor_flow"
    if haskey(data, "initial_compressor_flow") && haskey(data, "compressor_flow")
        (throw(ICException("both keys compressor_flow  and initial_compressor_flow present")))
    end

    compressor_flow_key = ""
    (haskey(data, "initial_compressor_flow")) && (compressor_flow_key = "initial_compressor_flow")
    (haskey(data, "compressor_flow")) && (compressor_flow_key = "compressor_flow")
    
    for (i, value) in get(data, compressor_flow_key, [])
        id = parse(Int64, i)
        ic[:compressor][id] = value
    end
    
    # initial condition for pipe flows (mandatory field)
    # initial conditions file should not have both "pipe_flow" and "initial_pipe_flow"
    if haskey(data, "initial_pipe_flow") && haskey(data, "pipe_flow")
        (throw(ICException("both keys pipe_flow  and initial_pipe_flow present")))
    end

    pipe_flow_key = haskey(data, "initial_pipe_flow") ? "initial_pipe_flow" : "pipe_flow"

    for (i, value) in get(data, pipe_flow_key, [])
        id = parse(Int64, i)
        L = data["pipes"][i]["length"]
        # if steady state (Number) else transient
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
    
    # if not given, computed later (this is not mandatory)
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
