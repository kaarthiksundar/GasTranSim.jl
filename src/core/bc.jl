function build_bc(data::Dict{String,Any})::Dict{Symbol,Any}
    bc = Dict{Symbol,Any}()

    bc[:node] = Dict() 
    bc[:compressor] = Dict()

    for (i, value) in get(data, "boundary_pslack", [])
        id = parse(Int64, i)
        t = value["time"]
        val = value["value"]
        bc[:node][id] = Dict( 
            "spl" => Spline1D(t, val, k=1), 
            "control_type" => pressure_control
        )
    end 

    for (i, value) in get(data, "boundary_nonslack_flow", [])
        id = parse(Int64, i)
        t = value["time"]
        val = value["value"] 
        bc[:node][id] = Dict(
            "spl" => Spline1D(t, val, k=1),
            "control_type" => flow_control
        )
    end 

    for (i, value) in get(data, "boundary_compressor", [])
        id = parse(Int64, i)
        t = value["time"]
        val = value["value"]
        bc_type = value["control_type"]
        bc[:compressor][id] = Dict(
            "spl" => CompressorControl(t, bc_type, val)
        )
    end 

    return bc
end 