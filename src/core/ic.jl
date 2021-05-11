function build_ic(data::Dict{String,Any})::Dict{Symbol,Any}
    ic = Dict{Symbol,Any}()

    ic[:node] = Dict() 
    ic[:pipe] = Dict("mass_flow" => Dict(), "pressure" => ())

    for (i, value) in get(data, "initial_nodal_pressure", [])
        id = parse(Int64, i)
        ic[:node][id] = value
    end 

    for (i, value) in get(data, "initial_pipe_flow", [])
        id = parse(Int64, i)
        L = data["pipes"][i]["length"]
        if (isa(value, Number))
            distance = [0.0, L]
            values = [value, value]
            ic[:pipe]["mass_flow"][id] = Spline1D(distance, values, k=1)
        else 
            distance = value["distance"]
            values = value["value"]
            ic[:pipe]["mass_flow"][id] = Spline1D(distance, values, k=1)
        end 
    end 

    for (i, value) in get(data, "initial_pipe_pressure", [])
        id = parse(Int64, i)
        distance = value["distance"]
        values = value["value"]
        ic[:pipe]["pressure"][id] = Spline1D(distance, values, k=1)
    end 
    return ic
end 

function add_initial_nodal_conditions_to_ref!(ref::Dict{Symbol,Any}, ic::Dict{Symbol,Any})

    (isempty(ic[:node])) && (@error "nodal pressure initial condition missing") 
    (isempty(ic[:pipe]["mass_flow"])) && (@error "pipe flow initial condition missing")

    for (i, value) in ic[:node]
        @assert i in keys(ref[:node])
        ref[:node][i]["initial_pressure"] = value
    end
    return
end
