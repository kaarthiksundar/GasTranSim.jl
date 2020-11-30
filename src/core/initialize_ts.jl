function initialize_simulator(file::AbstractString)::TransientSimulator
    data = parse_json(file)
    return initialize_simulator(data)
end 

function initialize_simulator(data::Dict{String,Any})::TransientSimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    ref = build_ref(data)
    bc = build_bc(data)
    pu_pressure_to_pu_density = x -> x 
    pu_density_to_pu_pressure = x -> x 
    pu_dp_drho = x-> x^0

    return TransientSimulator(data, 
        ref, 
        Dict{String,Any}(), 
        nominal_values, 
        params, 
        bc, 
        pu_pressure_to_pu_density, 
        pu_density_to_pu_pressure,
        pu_dp_drho) 
end 