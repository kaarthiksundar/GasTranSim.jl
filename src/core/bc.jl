function build_bc(data::Dict{String,Any})::Dict{Symbol,Any}
    bc = Dict{Symbol,Any}()

    bc[:node] = Dict() 
    bc[:compressor] = Dict()

    for (i, value) in get(data, "boundary_pslack", [])
        id = parse(Int64, i)
        time = value["time"]
        values = value["value"]
        bc[:node][id] = Dict( 
            "spl" => Spline1D(time, values, k=1), 
            "control_type" => pressure
        )
    end 

    for (i, value) in get(data, "boundary_nonslack_flow", [])
        id = parse(Int64, i)
        time = value["time"]
        values = value["value"]
        bc[:node][id] = Dict(
            "spl" => Spline1D(time, values, k=1),
            "control_type" => flow
        )
    end 

    for (i, value) in get(data, "boundary_compressor", [])
        id = parse(Int64, i)
        time = value["time"]
        values = value["value"]
        bc_type = value["control_type"]
        bc[:compressor][id] = Dict(
            "spl" => CompressorControl(time, bc_type, values)
        )
    end 

    return bc
end 