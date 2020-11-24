function si_to_english!(data::Dict{String,Any},
    params::Dict{Symbol,Any}, nominal_values::Dict{Symbol,Any})

    rescale_mass_flow = x -> x / get_mmscfd_to_kgps_conversion_factor(params)
    rescale_mass_flux = x -> x / get_mmscfd_to_kgps_conversion_factor(params)
    rescale_time = x -> x 
    rescale_pressure = x -> pascal_to_psi(x)
    rescale_length = x -> m_to_miles(x)
    rescale_density = x -> x
    rescale_diameter = x -> m_to_inches(x)

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