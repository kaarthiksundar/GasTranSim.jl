function get_eos(eos::Symbol)
    (eos == :ideal) &&
        (return _ideal_coeffs, _pressure_to_density_ideal, _density_to_pressure_ideal)
    (eos == :simple_cnga) && (
        return _simple_cnga_coeffs,
        _pressure_to_density_simple_cnga,
        _density_to_pressure_simple_cnga
    )
    (eos == :full_cnga) && (
        return _full_cnga_coeffs,
        _pressure_to_density_full_cnga,
        _density_to_pressure_full_cnga
    )

end

function _ideal_coeffs(
    nominal_values::Dict{Symbol,Any},
    params::Dict{Symbol,Any},
)::Tuple{Float64,Float64}
    return nominal_values[:euler_num], 0.0 # dimensionless
end

function _simple_cnga_coeffs(
    nominal_values::Dict{Symbol,Any},
    params::Dict{Symbol,Any},
)::Tuple{Float64,Float64}
    p0 = nominal_values[:pressure]
    Euler_num = nominal_values[:euler_num]
    b1 = 1.00300865  # dimensionless
    b2 = 2.96848838e-8 # units 1/pressure
    return Euler_num * b1, Euler_num * p0 * b2
end

function _full_cnga_coeffs(
    nominal_values::Dict{Symbol,Any},
    params::Dict{Symbol,Any},
)::Tuple{Float64,Float64}
    p0, p_atm = nominal_values[:pressure], 101350.0
    Euler_num = nominal_values[:euler_num]
    G, T = params[:gas_specific_gravity], params[:temperature]
    a1, a2, a3 = 344400.0, 1.785, 3.825 # all are dimensionless

    b1 = 1.0 + a1 * (p_atm/6894.75729) * (10 ^ (a2 * G)) / (1.8 * T) ^ a3 # dimensionless
    b2 = a1 * (10.0 ^ (a2 * G)) / (6894.75729 * (1.8 * T)^a3) # units 1/pressure
    return Euler_num * b1, Euler_num * p0 * b2
end

function _pressure_to_density_ideal(
    p,
    nominal_values::Dict{Symbol,Any},
    params::Dict{Symbol,Any},
)
    Euler_num = nominal_values[:euler_num]
    return Euler_num * p
end

function _density_to_pressure_ideal(
    rho,
    nominal_values::Dict{Symbol,Any},
    params::Dict{Symbol,Any},
)
    Euler_num = nominal_values[:euler_num]
    return rho/Euler_num
end

function _pressure_to_density_simple_cnga(
    p,
    nominal_values::Dict{Symbol,Any},
    params::Dict{Symbol,Any},
)
    b1, b2 = _simple_cnga_coeffs(nominal_values, params)
    return (p .* (b1 .+ b2 * p))
end

function _density_to_pressure_simple_cnga(
    rho,
    nominal_values::Dict{Symbol,Any},
    params::Dict{Symbol,Any},
)
    b1, b2 = _simple_cnga_coeffs(nominal_values, params)
    return (-b1 .+ sqrt.(b1 * b1 .+ 4 * b2 * rho)) / (2.0 * b2)
end

function _pressure_to_density_full_cnga(
    p,
    nominal_values::Dict{Symbol,Any},
    params::Dict{Symbol,Any},
)
    b1, b2 = _full_cnga_coeffs(nominal_values, params)
    return (p .* (b1 .+ b2 * p))
end

function _density_to_pressure_full_cnga(
    rho,
    nominal_values::Dict{Symbol,Any},
    params::Dict{Symbol,Any},
)
    b1, b2 = _full_cnga_coeffs(nominal_values, params)
    return (-b1 .+ sqrt.(b1 * b1 .+ 4 * b2 * rho)) / (2.0 * b2)
end
