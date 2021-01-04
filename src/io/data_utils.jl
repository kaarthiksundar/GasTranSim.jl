function process_data!(data::Dict{String,Any})
    nominal_values = Dict{Symbol,Any}()
    params = Dict{Symbol,Any}()

    params_exhaustive = ["temperature", "gas_specific_gravity",
        "specific_heat_capacity_ratio", "units", "t_0", "t_f", "dt",
        "courant_number", "nodal_output_dt", "pipe_output_flag",
        "pipe_output_dt", "pipe_output_dx"
    ]
    defaults_exhaustive = [288.706, 0.6, 1.4, 0, 0.0, 3600.0, 1.0, 0.95, 600.0,
        0, 3600.0, 1000.0
    ]
    input_params = data["input_param"]
    key_map = Dict{String,String}()
    for k in keys(input_params)
        occursin("Temperature", k) && (key_map["temperature"] = k)
        occursin("Gas", k) && (key_map["gas_specific_gravity"] = k)
        occursin("Specific heat", k) &&
            (key_map["specific_heat_capacity_ratio"] = k)
        occursin("units", k) && (key_map["units"] = k)
        occursin("Initial time", k) && (key_map["t_0"] = k)
        occursin("Final time", k) && (key_map["t_f"] = k)
        occursin("Discretization", k) && (key_map["dt"] = k)
        occursin("Courant", k) && (key_map["courant_number"] = k)
        occursin("Output nodal dt", k) && (key_map["nodal_output_dt"] = k)
        occursin("Output pipe data", k) && (key_map["pipe_output_flag"] = k)
        occursin("Output pipe dt", k) && (key_map["pipe_output_dt"] = k)
        occursin("dx", k) && (key_map["pipe_output_dx"] = k)

    end

    # add "area" key to pipes in data
    for (i, pipe) in get(data, "pipes", [])
        pipe["area"] = pi * pipe["diameter"] * pipe["diameter"] * 0.25
    end

    # populating parameters
    for i in 1:length(params_exhaustive)
        param = params_exhaustive[i]
        default = defaults_exhaustive[i]
        if param == "units"
            if haskey(key_map, param)
                value = Int(input_params[key_map[param]])
                if (value == 0)
                    params[:units] = 0
                    params[:is_si_units] = 1
                    params[:is_english_units] = 0
                    params[:is_per_unit] = 0
                else
                    params[:units] = 0
                    params[:is_is_units] = 0
                    params[:is_english_units] = 1
                    params[:is_per_unit] = 0
                end
            else
                params[:units] = 0
                params[:is_si_units] = 1
                params[:is_english_units] = 0
                params[:is_per_unit] = 0
            end
            continue
        end
        key = get(key_map, param, false)
        if key != false
            value = input_params[key]
            params[Symbol(param)] = value
        else
            params[Symbol(param)] = default
        end
    end

    # other parameter calculations
    # universal gas constant (J/mol/K)
    params[:R] = 8.314
    # molecular mass of natural gas (kg/mol): M_g = M_a * G
    params[:gas_molar_mass] = 0.02896 * params[:gas_specific_gravity]
    # sound speed (m/s): v = sqrt(R_g * T); R_g = R/M_g = R/M_a/G; R_g is specific gas constant; g-gas, a-air
    params[:sound_speed] = sqrt(params[:R] * params[:temperature] / params[:gas_molar_mass])

    nominal_values[:length] = 5000.0
    nominal_values[:area] = 1.0
    nominal_values[:pressure] = 3500000.0 # 507.63 psi
    nominal_values[:density] = nominal_values[:pressure] / params[:sound_speed] / params[:sound_speed]
    nominal_values[:mass_flux] = nominal_values[:density] * params[:sound_speed]
    nominal_values[:time] = nominal_values[:length] / params[:sound_speed]
    nominal_values[:mass_flow] = nominal_values[:mass_flux] * nominal_values[:area]

    return params, nominal_values
end
