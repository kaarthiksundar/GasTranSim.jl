function initialize_simulator(file::AbstractString)::TransientSimulator
    data = parse_json(file)
    return initialize_simulator(data)
end 

function initialize_simulator(data::Dict{String,Any})::TransientSimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    ref = build_ref(data)
    bc = build_bc(data)

    return TransientSimulator(data, 
        ref, 
        Dict{String,Any}(), 
        nominal_values, 
        params, 
        bc) 
    
end 