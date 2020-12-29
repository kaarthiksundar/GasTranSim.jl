function pu_to_si!(data::Dict{String,Any},
    params::Dict{Symbol,Any}, nominal_values::Dict{Symbol,Any})

    rescale_mass_flow = x -> x * nominal_values[:mass_flow]
    rescale_mass_flux = x -> x * nominal_values[:mass_flux]
    rescale_time = x -> x * nominal_values[:time]
    rescale_pressure = x -> x * nominal_values[:pressure]
    rescale_length = x -> x * nominal_values[:length]
    rescale_density = x -> x * nominal_values[:density]
    rescale_diameter = x -> x * nominal_values[:length]

    function rescale_compressor_boundary_conditions!(type, value)
        @assert length(type) == length(value)
        for i in 1:length(type)
            (type[i] == 1) && (value[i] = rescale_pressure(value[i]))
            (type[i] == 2) && (value[i] = rescale_mass_flow(value[i]))
        end 
    end 

    rescale_functions = [rescale_mass_flow, rescale_mass_flux, 
        rescale_time, rescale_pressure, rescale_length, rescale_density, 
        rescale_diameter, rescale_compressor_boundary_conditions!]
    
    _rescale_data!(data, params, rescale_functions)
end 

function english_to_si!(data::Dict{String,Any},
    params::Dict{Symbol,Any}, nominal_values::Dict{Symbol,Any})

    rescale_mass_flow = x -> x * get_mmscfd_to_kgps_conversion_factor(params)
    rescale_mass_flux = x -> x * get_mmscfd_to_kgps_conversion_factor(params)
    rescale_time = x -> x 
    rescale_pressure = x -> psi_to_pascal(x)
    rescale_length = x -> miles_to_m(x)
    rescale_density = x -> x
    rescale_diameter = x -> inches_to_m(x)

    function rescale_compressor_boundary_conditions!(type, value)
        @assert length(type) == length(value)
        for i in 1:length(type)
            (type[i] == 1) && (value[i] = rescale_pressure(value[i]))
            (type[i] == 2) && (value[i] = rescale_mass_flow(value[i]))
        end 
    end 

    rescale_functions = [rescale_mass_flow, rescale_mass_flux, 
        rescale_time, rescale_pressure, rescale_length, rescale_density, 
        rescale_diameter, rescale_compressor_boundary_conditions!]
    
    _rescale_data!(data, params, rescale_functions)
end 

