function get_eos(nominal_values::Dict{Symbol,Any}, params::Dict{Symbol,Any}, eos::Symbol)
    (eos == :ideal) &&
        (return _pressure_to_density_ideal, _density_to_pressure_ideal)
    (eos == :simple_cnga) &&
        (return _pressure_to_density_simple_cnga, _density_to_pressure_simple_cnga)
    (eos == :full_cnga) &&
        (return _pressure_to_density_full_cnga, _density_to_pressure_full_cnga)

end

function _pressure_to_density_ideal(p, nominal_values::Dict{Symbol,Any}, params::Dict{Symbol,Any})
    return p  # since non-dimensional p = rho
end

function _density_to_pressure_ideal(rho, nominal_values::Dict{Symbol,Any}, params::Dict{Symbol,Any})
    return rho  # since non-dimensional p = rho
end

function _pressure_to_density_simple_cnga(p, nominal_values::Dict{Symbol,Any}, params::Dict{Symbol,Any})
    p0, rho0 = nominal_values[:pressure], nominal_values[:density]
    b1, b2 = 1.00300865, 2.96848838e-8
    RgT = nominal_values[:sound_speed]*nominal_values[:sound_speed]
    return (p0 * p .* (b1 .+ b2 * p0 * p ) / RgT ) / rho0
end

function _density_to_pressure_simple_cnga(rho, nominal_values::Dict{Symbol,Any}, params::Dict{Symbol,Any})
    p0, rho0 = nominal_values[:pressure], nominal_values[:density]
    b1, b2 = 1.00300865, 2.96848838e-8
    RgT = nominal_values[:sound_speed]*nominal_values[:sound_speed]

    return (-b1 .+ sqrt.(b1 * b1 .+ 4 * b2 * RgT * rho0 * rho) ) / (2.0 * b2) / p0
end

function _pressure_to_density_full_cnga(p, nominal_values::Dict{Symbol,Any}, params::Dict{Symbol,Any})
    p0, rho0, p_atm = nominal_values[:pressure], nominal_values[:density], 101350.0
    G, T = params[:gas_specific_gravity], params[:temperature]
    a1, a2, a3 = 344400.0, 1.785, 3.825
    b1 = 1.0 + a1 * (p_atm/6894.75729) * ( 10 ^ (a2 * G) ) / (1.8 * T) ^ a3
    b2 = a1 * (10.0 ^ (a2 * G) ) / (6894.75729 * (1.8 * T)^a3 )
    RgT = nominal_values[:sound_speed]*nominal_values[:sound_speed]

    return (p0 * p .* (b1 .+ b2 * p0 * p ) / (RgT) ) / rho0
end

function _density_to_pressure_full_cnga(rho, nominal_values::Dict{Symbol,Any}, params::Dict{Symbol,Any})
    p0, rho0, p_atm = nominal_values[:pressure], nominal_values[:density], 101350
    G, T = params[:gas_specific_gravity], params[:temperature]
    a1, a2, a3 = 344400, 1.785, 3.825
    b1 = 1 + a1 * (p_atm/6894.75729) * ( 10 ^ (a2 * G) ) / (1.8 * T) ^ a3
    b2 = a1 * (10 ^ (a2 * G) ) / (6894.75729 * (1.8 * T)^a3 )
    RgT = nominal_values[:sound_speed]*nominal_values[:sound_speed]

    return (-b1 .+ sqrt.(b1 * b1 .+ 4 * b2 * RgT * rho0 * rho) ) / (2 * b2) / p0
end
