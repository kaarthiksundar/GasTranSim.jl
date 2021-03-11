
@inline psi_to_pascal(psi) = psi * 6894.75729
@inline pascal_to_psi(pascal) = pascal / 6894.75729
@inline km_to_m(km) = km * 1000.0
@inline m_to_km(m) = m / 1000.0
@inline hp_to_watts(hp) = hp * 745.7
@inline watts_to_hp(watts) = watts / 745.7
@inline miles_to_m(miles) = miles * 1609.64
@inline m_to_miles(m) = m / 1609.64
@inline inches_to_m(inches) = inches * 0.0254
@inline m_to_inches(m) = m / 0.0254
@inline sq_inches_to_sq_m(sq_inches) = sq_inches * 0.0254 * 0.0254
@inline sq_m_to_sq_inches(sq_m) = sq_m / (0.0254 * 0.0254)
@inline get_universal_R() = 8.314
@inline get_universal_R(params::Dict{Symbol,Any}) = get(params, :R, get_universal_R())
@inline get_gas_specific_gravity(params::Dict{Symbol,Any}) = get(params, :gas_specific_gravity, 0.6)
@inline get_temperature(params::Dict{Symbol,Any}) = get(params, :temperature, 288.7060)
@inline get_sound_speed(params::Dict{Symbol,Any}) = get(params, :sound_speed, 371.6643)
@inline get_molecular_mass_of_air() = 0.02896
@inline get_one_atm_in_pascal() = 101325.0
@inline get_one_atm_in_psi() = pascal_to_psi(101325)

"""
Convering mmsfcd to kgps (standard volumetric flow rate to mass flow rate)
"""

function get_mmscfd_to_kgps_conversion_factor(params::Dict{Symbol,Any})::Number
    standard_pressure = get_one_atm_in_pascal()
    R = get_universal_R(params)
    standard_temperature = get_temperature(params)
    cubic_ft_to_cubic_m = 0.02832
    volumetric_flow_rate_in_si = cubic_ft_to_cubic_m * 1e6 / 86400.0
    molecular_mass_of_gas = get_gas_specific_gravity(params) * get_molecular_mass_of_air()
    density_at_standard_conditions = standard_pressure * molecular_mass_of_gas / standard_temperature / R
    return density_at_standard_conditions * volumetric_flow_rate_in_si
end

get_kgps_to_mmscfd_conversion_factor(params::Dict{Symbol,Any})::Number = 1 / get_mmscfd_to_kgps_conversion_factor(params)

function get_data_units(rescale_functions)::Dict{Symbol,Any}

    units = Dict{Symbol,Any}()
    rescale_mass_flow = rescale_functions[1]
    rescale_mass_flux = rescale_functions[2]
    rescale_time = rescale_functions[3]
    rescale_pressure = rescale_functions[4]
    rescale_length = rescale_functions[5]
    rescale_density = rescale_functions[6]
    rescale_diameter = rescale_functions[7]
    rescale_area = rescale_functions[8]
    params_units = Dict{String,Function}(
        "t_0" => rescale_time,
        "t_f" => rescale_time, 
        "dt" => rescale_time, 
        "output_dx" => rescale_length,
        "output_dt" => rescale_time
    )

    node_units = Dict{String,Function}(
        "min_pressure" => rescale_pressure,
        "max_pressure" => rescale_pressure, 
        "min_injection" => rescale_mass_flow, 
        "max_injection" => rescale_mass_flow, 
    )

    pipe_units = Dict{String,Function}(
        "diameter" => rescale_diameter,
        "length" => rescale_length,
        "area" => rescale_area,
    )

    compressor_units = Dict{String,Function}(
        "max_power" => rescale_mass_flow, 
        "min_flow" => rescale_mass_flow, 
        "max_flow" => rescale_mass_flow, 
    )

    initial_flow_units = Dict{String,Function}(
        "initial_nodal_flow" => rescale_mass_flow,
        "initial_pipe_flow" => rescale_mass_flow, 
        "initial_compressor_flow" => rescale_mass_flow
    )

    initial_pressure_units = Dict{String,Function}(
        "initial_pipe_pressure_in" => rescale_pressure, 
        "initial_pipe_pressure_out" => rescale_pressure, 
        "initial_compressor_pressure_in" => rescale_pressure, 
        "initial_compressor_pressure_out" => rescale_pressure, 
        "initial_nodal_pressure" => rescale_pressure
    )

    boundary_flow_units = Dict{String,Function}(
        "boundary_nonslack_flow" => rescale_mass_flow 
    )

    boundary_pressure_units = Dict{String,Function}(
        "boundary_pslack" => rescale_pressure
    )

    node_disruption_units = Dict{String,Function}(
        "disrupt_time" => rescale_time,
        "valve_delay" => rescale_time,
    )

    pipe_disruption_units = Dict{String,Function}(
        "disrupt_time" => rescale_time,
        "valve_delay" => rescale_time,
        "distance" => rescale_length
    )

    units[:params_units] = params_units
    units[:node_units] = node_units
    units[:pipe_units] = pipe_units
    units[:compressor_units] = compressor_units
    units[:initial_flow_units] = initial_flow_units
    units[:initial_pressure_units] = initial_pressure_units
    units[:boundary_flow_units] = boundary_flow_units 
    units[:boundary_pressure_units] = boundary_pressure_units
    units[:node_disruption_units] = node_disruption_units 
    units[:pipe_disruption_units] = pipe_disruption_units

    return units
end 

function _rescale_data!(data::Dict{String,Any}, 
    params::Dict{Symbol,Any}, rescale_functions::Vector{Function})

    units = get_data_units(rescale_functions)
    params_units = units[:params_units]
    node_units = units[:node_units]
    pipe_units = units[:pipe_units]
    compressor_units = units[:compressor_units] 
    initial_flow_units = units[:initial_flow_units]
    initial_pressure_units = units[:initial_pressure_units]
    boundary_flow_units = units[:boundary_flow_units]
    boundary_pressure_units = units[:boundary_pressure_units]
    node_disruption_units = units[:node_disruption_units]
    pipe_disruption_units = units[:pipe_disruption_units]

    rescale_mass_flow = rescale_functions[1]
    rescale_mass_flux = rescale_functions[2]
    rescale_time = rescale_functions[3]
    rescale_pressure = rescale_functions[4]
    rescale_length = rescale_functions[5]
    rescale_density = rescale_functions[6]
    rescale_diameter = rescale_functions[7]
    rescale_area = rescale_functions[8]
    rescale_boundary_compressor = rescale_functions[9]

    for (param, f) in params_units         
        value = params[Symbol(param)]
        params[Symbol(param)] = f(value)
    end 

    for (i, node) in get(data, "nodes", [])
        for (param, f) in node_units
            (!haskey(node, param)) && (continue)
            value = node[param]
            node[param] = f(value)
        end 
    end 

    for (i, pipe) in get(data, "pipes", [])
        for (param, f) in pipe_units
            (!haskey(pipe, param)) && (continue)
            value = pipe[param]
            pipe[param] = f(value)
        end 
    end 
    
    for (i, compressor) in get(data, "compressors", [])
        for (param, f) in compressor_units
            (!haskey(Dict(compressor), param)) && (continue)
            value = compressor[param]
            compressor[param] = f(value)
        end 
    end 

    for (param, f) in merge(initial_flow_units, initial_pressure_units)
        for (i, value) in get(data, param, [])
            data[param][i] = f(value)
        end 
    end 

    for (param, f) in merge(boundary_flow_units, boundary_pressure_units)
        for (i, value) in get(data, param, [])
            data[param][i]["time"] = rescale_time.(value["time"])
            data[param][i]["value"] = f.(value["value"])
        end 
    end 

    for (i, value) in get(data, "boundary_compressor", []) 
        data["boundary_compressor"][i]["time"] = rescale_time.(value["time"])
        rescale_boundary_compressor(value["control_type"], value["value"])
    end 

    for (i, disruption) in get(data, "disruption", [])
        if (i == "node_disruption")
            for (param, f) in node_disruption_units
                for (i, node_disruption) in disruption
                    value = node_disruption[param]
                    node_disruption[param] = f(value)
                end 
            end 
        end 
        if (i == "pipe_disruption")
            for (param, f) in pipe_disruption_units
                for (i, pipe_disruption) in disruption
                    value = pipe_disruption[param]
                    pipe_disruption[param] = f(value)
                end 
            end 
        end 
    end 
end 