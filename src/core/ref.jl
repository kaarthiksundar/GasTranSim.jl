function build_ref(data::Dict{String,Any})::Dict{Symbol,Any}
    component_types = ["nodes", "pipes", "compressors", "gnodes"]
    initial_conditions = Dict(
        "nodes" => ["initial_nodal_flow", "initial_nodal_pressure"],
        "pipes" => ["initial_pipe_pressure_in", "initial_pipe_pressure_out", "initial_pipe_flow"],
        "compressors" => ["initial_compressor_flow", "initial_compressor_pressure_in", 
            "initial_compressor_pressure_out", "initial_compressor_ratio"],
        "gnodes" => []
    )
    boundary_conditions = Dict(
        "nodes" => ["boundary_nonslack_flows", "boundary_pslack"],
        "pipes" => [],
        "compressors" => ["boundary_compressor"], 
        "gnodes" => []
    )

end 