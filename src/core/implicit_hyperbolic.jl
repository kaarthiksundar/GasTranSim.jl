
function implicit_hyperbolic_step!(ts::TransientSimulator, run_type::Symbol)
    advance_junction_pressures!(ts, :implicit_hyperbolic, run_type)
    return
end


function initialize_pipe_grid!(ts::TransientSimulator, ::Val{:implicit_hyperbolic})
    for (key, pipe) in ref(ts, :pipe)
        num_cvs = get(params(ts), :pipe_segments, 20)
        ref(ts, :pipe, key)["num_discretization_points"] = num_cvs
        ref(ts, :pipe, key)["dx"] = pipe["length"] / (num_cvs + 1)
        ref(ts, :pipe, key)["density_profile"] = zeros(Float64, num_cvs)
        ref(ts, :pipe, key)["mass_flux_profile"] = zeros(Float64,num_cvs)

    end
    return
end

function initialize_pipe_state!(ts::TransientSimulator, ::Val{:implicit_hyperbolic})
    is_steady = false
    # We are assuming initial pipe pressures will not be provided only for steady
    # initial conditions. In the unsteady case we assume initial pipe pressures are given.
    if isempty(ts.initial_conditions[:pipe]["pressure"])
        @info "Pipes do not have initial pressure profile, will be computed assuming steady state flow"
        is_steady = true
    end
    for (key, pipe) in ref(ts, :pipe)
        area = pipe["area"]
        n = pipe["num_discretization_points"]
        dx = pipe["dx"]
        L = pipe["length"]
        fr_node = pipe["fr_node"]
        to_node = pipe["to_node"]
        x_centers = ((collect(1:n) .- 0.5) .* dx)
        if is_steady
            initial_mass_flux = initial_pipe_mass_flow(ts, key)(0.0) / area
            fill!(pipe["mass_flux_profile"], initial_mass_flux)
            initial_fr_pressure = ref(ts, :node, fr_node, "pressure")
            initial_to_pressure = ref(ts, :node, to_node, "pressure")
            density_at_first_sq = get_density(ts, initial_fr_pressure)^2
            density_at_last_sq = get_density(ts, initial_to_pressure)^2
            pipe["density_profile"][1:n] = [
                sqrt(density_at_first_sq + (density_at_last_sq - density_at_first_sq) * (x / L))
                for x in x_centers
            ]
            
            pipe["fr_minus_mass_flux"] = initial_mass_flux # first cell center (x = dx/2)
            pipe["to_minus_mass_flux"] = initial_mass_flux # last cell center (x = L - dx/2)
            pipe["fr_mass_flux"] = initial_mass_flux # 0
            pipe["to_mass_flux"] = initial_mass_flux # L
        else
            flow_spl = initial_pipe_mass_flow(ts, key)
            pressure_spl = initial_pipe_pressure(ts, key)
            pipe["mass_flux_profile"][1:n] = [flow_spl(x) for x in x_centers] ./ area
            pipe["fr_minus_mass_flux"] = pipe["mass_flux_profile"][1]
            pipe["to_minus_mass_flux"] = pipe["mass_flux_profile"][end]
            pipe["fr_mass_flux"] = flow_spl(0.0) / area
            pipe["to_mass_flux"] = flow_spl(L) / area
            pipe["density_profile"][1:n] = [get_density(ts, pressure_spl(x)) for x in x_centers]
        end
        pipe["phi"] = zeros(n)
        pipe["rho"] = zeros(n)
    end
    return
end



function _solve_pipe_state_hyperbolic!(
    ts::TransientSimulator,
    pipe_id::Int64,
    rho_from::T,
    rho_to::T)::Vector{T} where {T<:Real}
    n = ref(ts, :pipe, pipe_id)["num_discretization_points"]
    x = zeros(T, 2 * n)
    x[1:n] = T.(ref(ts, :pipe, pipe_id)["density_profile"])
    x[n+1:2*n] = T.(ref(ts, :pipe, pipe_id)["mass_flux_profile"])

    rho_old = T.(ref(ts, :pipe, pipe_id)["density_profile"])
    phi_old = T.(ref(ts, :pipe, pipe_id)["mass_flux_profile"])
    area = T(ref(ts, :pipe, pipe_id)["area"])


    residual_fun! = (r, x) -> _pipe_residual_hyperbolic!(
        r,
        ts,
        pipe_id,
        x,
        rho_old,
        phi_old,
        rho_from,
        rho_to
    )
    Jacobian_fun! = (J, x) -> _pipe_jacobian_hyperbolic!(
        J,
        ts,
        pipe_id,
        x,
        rho_from,
        rho_to
    )


    x, converged, _, res_norm = solve_newton_basic!(x, residual_fun!, Jacobian_fun!, tol = 1e-6, max_iter = 100)

    
    converged || throw(DomainError(res_norm, "Newton solver did not converge for pipe $pipe_id"))

    # phi_from, phi_to = get_bry_flux_char_extrapolation(ts, rho_from, rho_to, x[1], x[n], x[n+1], x[n+n])
    phi_from, phi_to = get_bdry_flux_2nd_order_char_extrapolation(ts, rho_from, rho_to, x[1], x[2], x[n-1], x[n], x[n+1], x[n+2], x[n+n-1], x[n+n])

    if T == Float64
        ref(ts, :pipe, pipe_id)["rho"] = x[1:n]
        ref(ts, :pipe, pipe_id)["phi"] = x[n+1:2*n]
        ref(ts, :pipe, pipe_id)["fr_mass_flux"] = phi_from
        ref(ts, :pipe, pipe_id)["to_mass_flux"] = phi_to
    end


    end_flows = T[area * phi_from, area * phi_to]
    
    return end_flows

end

function get_bry_flux_char_extrapolation(ts::TransientSimulator, rho_from::T, rho_to::T, rho_1::T, rho_n::T, phi_1::T, phi_n::T)::Tuple{T, T} where {T<:Real}

    
    pot_from = characteristic_potential(ts, rho_from)
    pot_1 = characteristic_potential(ts, rho_1)
    pot_to = characteristic_potential(ts, rho_to)
    pot_n = characteristic_potential(ts, rho_n)

    if params(ts, :inertial_flag) == false
        phi_from = phi_1  - pot_1 + pot_from
        phi_to = phi_n   + pot_n - pot_to
    else
        phi_from = rho_from * (phi_1 / rho_1  - pot_1 + pot_from)
        phi_to = rho_to * (phi_n / rho_n   + pot_n - pot_to)
    end

    return phi_from, phi_to
end

function get_bry_flux_matrices_for_char_extrapolation(ts::TransientSimulator, rho_from::T, rho_to::T, rho_1::T, rho_n::T,phi_1::T, phi_n::T)::Tuple{Matrix{T}, Matrix{T}} where {T<:Real}

    
    if params(ts, :inertial_flag) == false
        Bdry_left_mat = [T(0) T(0); -characteristic_potential_prime(ts, rho_1)  T(1)]
        Bdry_right_mat = [T(0) T(0);  characteristic_potential_prime(ts, rho_n)  T(1)]
    else
        # for Jacobian
        Bdry_left_mat = [T(0) T(0); rho_from * ( -phi_1 / (rho_1 ^ 2) - characteristic_potential_prime(ts, rho_1))  rho_from * (1.0 / rho_1)]
        Bdry_right_mat = [T(0) T(0); rho_to * ( -phi_n / (rho_n ^ 2) + characteristic_potential_prime(ts, rho_n) ) rho_to * ( 1.0 / rho_n)]
    end

    return Bdry_left_mat, Bdry_right_mat
end

function get_bdry_flux_2nd_order_char_extrapolation(ts::TransientSimulator, rho_from::T, rho_to::T, rho_1::T, rho_2::T, rho_n_minus_1::T, rho_n::T, phi_1::T, phi_2::T, phi_n_minus_1::T, phi_n::T)::Tuple{T, T} where {T<:Real}

    pot_from = characteristic_potential(ts, rho_from)
    pot_1 = characteristic_potential(ts, rho_1)
    pot_2 = characteristic_potential(ts, rho_2)

    pot_to = characteristic_potential(ts, rho_to)
    pot_n = characteristic_potential(ts, rho_n)
    pot_n_minus_1 = characteristic_potential(ts, rho_n_minus_1)

    

    if params(ts, :inertial_flag) == false
        phi_from = (3/2) * (phi_1   - pot_1) - (1/2) * (phi_2  - pot_2) + pot_from 
        phi_to = (3/2) * (phi_n  + pot_n) - (1/2) * (phi_n_minus_1  + pot_n_minus_1) - pot_to
    else 
        phi_from = rho_from * ( (3/2) * (phi_1 / rho_1   - pot_1) - (1/2) * (phi_2 / rho_2 - pot_2) + pot_from )
        phi_to = rho_to * ( (3/2) * (phi_n / rho_n + pot_n) - (1/2) * (phi_n_minus_1 / rho_n_minus_1  + pot_n_minus_1) - pot_to ) 
    end
    

    return phi_from, phi_to
end

function get_bdry_flux_matrices_for_2nd_order_char_extrapolation(ts::TransientSimulator, rho_from::T, rho_to::T, rho_1::T, rho_2::T, rho_n_minus_1::T, rho_n::T, phi_1::T, phi_2::T, phi_n_minus_1::T, phi_n::T)::Tuple{Matrix{T}, Matrix{T}, Matrix{T}, Matrix{T}} where {T<:Real}

    Bdry_left_mat_1 = zeros(T, 2, 2)
    Bdry_left_mat_2 = zeros(T, 2, 2)
    Bdry_right_mat_n = zeros(T, 2, 2)
    Bdry_right_mat_n_minus_1 = zeros(T, 2, 2)

    if params(ts, :inertial_flag) == false
        Bdry_left_mat_1 = [
            T(0) T(0);
            -T(3 / 2) * characteristic_potential_prime(ts, rho_1) T(3 / 2)
        ]
        Bdry_left_mat_2 = [
            T(0) T(0);
            T(1 / 2) * characteristic_potential_prime(ts, rho_2) -T(1 / 2)
        ]
        Bdry_right_mat_n_minus_1 = [
            T(0) T(0);
            -T(1 / 2) * characteristic_potential_prime(ts, rho_n_minus_1) -T(1 / 2)
        ]
        Bdry_right_mat_n = [
            T(0) T(0);
            T(3 / 2) * characteristic_potential_prime(ts, rho_n) T(3 / 2)
        ]
    else

        Bdry_left_mat_1[2, 1] =  rho_from * T(3/2) * ( -phi_1 / ( rho_1 ^ 2 ) - characteristic_potential_prime(ts, rho_1) )
        Bdry_left_mat_1[2, 2] = rho_from * T(3/2) * ( 1.0/ rho_1)


        Bdry_left_mat_2[2, 1] = -rho_from * T(1/2) * ( -phi_2 / (rho_2 ^ 2 )  - characteristic_potential_prime(ts, rho_2) )
        Bdry_left_mat_2[2, 2] = -rho_from * T(1/2) * ( 1.0 / rho_2 )

        Bdry_right_mat_n_minus_1[2, 1] = -rho_to * T(1/2) * ( -phi_n_minus_1 / (rho_n_minus_1 ^ 2 )  + characteristic_potential_prime(ts, rho_n_minus_1) )
        Bdry_right_mat_n_minus_1[2, 2] =  -rho_to * T(1/2) *  ( 1.0 /rho_n_minus_1)

        Bdry_right_mat_n[2, 1] = rho_to * T(3/2) * ( -phi_n / (rho_n ^ 2 ) + characteristic_potential_prime(ts, rho_n) )
        Bdry_right_mat_n[2, 2] = rho_to * T(3/2) * (  1.0 / rho_n )
    end

    return Bdry_left_mat_1, Bdry_left_mat_2, Bdry_right_mat_n_minus_1, Bdry_right_mat_n
end


"""
Residual for implicit hyperbolic pipe solve, matching solve_pipe_state! layout.

Unknown/state ordering:
- x[1:n]     -> rho
- x[n+1:2n]  -> phi

Residual ordering:
- r[1:n]     -> mass equations
- r[n+1:2n]  -> momentum equations
"""
function _pipe_residual_hyperbolic!(
    r::Vector{T},
    ts::TransientSimulator,
    pipe_id::Int64,
    x::Vector{T},
    rho_old::AbstractVector{T},
    phi_old::AbstractVector{T},
    rho_from::T,
    rho_to::T) where {T<:Real}

    n = div(length(x) , 2)

    rho = x[1:n]
    phi = x[n+1:2*n]
    @assert length(phi) == n
    @assert length(rho_old) == n
    @assert length(r) == 2n
    @assert n > 2

    fill!(r, zero(T))

    dx = ref(ts, :pipe, pipe_id)["dx"]
    dt = params(ts, :dt)
    mu =  dt / dx
    drag_coeff = ref(ts, :pipe, pipe_id, "friction_factor") /
        (2 * ref(ts, :pipe, pipe_id, "diameter"))
    g_coeff = ref(ts, :pipe, pipe_id, "sin_incline")  / (nominal_values(ts, :froude_num))^2
    p_coeff = nominal_values(ts, :euler_num) / (nominal_values(ts, :mach_num))^2

    r_local = zeros(T, 2)

    

    phi_from, phi_to = get_bdry_flux_2nd_order_char_extrapolation(ts, rho_from, rho_to, rho[1], rho[2], rho[n-1], rho[n], phi[1], phi[2], phi[n-1], phi[n])


    for i = 1 : n
        fill!(r_local, zero(T))

        U =  [ rho[i], phi[i] ]
        U_old = [rho_old[i], phi_old[i]]
        r_local += U - U_old 
        r_local += - dt * source(U, drag_coeff, g_coeff)
        if i == 1
            flux_i_plus_half  = numerical_flux(ts, U, [rho[i+1], phi[i+1]], p_coeff) 
            flux_i_minus_half = numerical_flux(ts, [rho_from, phi_from], U, p_coeff) 
        elseif i == n
            flux_i_plus_half  = numerical_flux(ts, U, [rho_to, phi_to], p_coeff) 
            flux_i_minus_half = numerical_flux(ts, [rho[i-1], phi[i-1]], U, p_coeff)    
        else
            flux_i_plus_half  = numerical_flux(ts, U, [ rho[i+1], phi[i+1] ], p_coeff) 
            flux_i_minus_half = numerical_flux(ts, [rho[i-1], phi[i-1]], U, p_coeff)
        end
        r_local +=  mu * ( flux_i_plus_half - flux_i_minus_half )
        assemble_local_residual!(ts, r, r_local, i, n)
    end
    
    return r
end


function rusanov_speed(ts::TransientSimulator, U_left::AbstractVector{T}, U_right::AbstractVector{T})::T where {T<:Real}
    epsilon = T(1e-6)
    vL = U_left[2] / U_left[1] # velocity on the left
    vR = U_right[2] / U_right[1] # velocity on the right
    cL = sqrt(get_pressure_prime(ts, U_left[1])) # speed of sound on the left
    cR = sqrt(get_pressure_prime(ts, U_right[1])) # speed of sound on the right
    if params(ts, :inertial_flag) == false
        sL =  cL
        sR =  cR
    else
        sL = smooth_abs(vL, epsilon) + cL
        sR = smooth_abs(vR, epsilon) + cR
    end

    return smooth_max(sL, sR, epsilon)
end

# function smooth_max(a::T, b::T, epsilon::T)::T where {T<:Real}
#     return 0.5 * (a + b + sqrt((a - b)^2 + epsilon^2))
# end

# function smooth_abs(x::T, epsilon::T)::T where {T<:Real}
#     return sqrt(x^2 + epsilon^2)
# end

function rusanov_speed_derivatives(ts::TransientSimulator, U_left::AbstractVector{T}, U_right::AbstractVector{T})::Tuple{Vector{T}, Vector{T}} where {T<:Real}
    
    epsilon = T(1e-6)
    vL = U_left[2] / U_left[1] # velocity on the left
    vR = U_right[2] / U_right[1] # velocity on the right
    cL = sqrt(get_pressure_prime(ts, U_left[1])) # speed of sound on the left
    cR = sqrt(get_pressure_prime(ts, U_right[1])) # speed of sound on the right

    if params(ts, :inertial_flag) == false
        sL = cL
        sR = cR
        t2_left = 0.0
        t2_right = 0.0
    else
        sL = smooth_abs(vL, epsilon) + cL
        sR = smooth_abs(vR, epsilon) + cR
        t2_left = vL / ( U_left[1] * smooth_abs(vL, epsilon) )
        t2_right = vR / ( U_right[1] * smooth_abs(vR, epsilon) )


    end

    wL  = T(1/2) * (1 + (sL - sR) / sqrt((sL - sR)^2 + epsilon^2) )
    wR  = T(1/2) * (1 - (sL - sR) / sqrt((sL - sR)^2 + epsilon^2) )

    ds_dU_left = zeros(T, 2)

    t1_left = get_pressure_double_prime(ts, U_left[1]) / (2 * cL)
    ds_dU_left[1] = wL * (t1_left - vL * t2_left)
    ds_dU_left[2] = wL * t2_left

    ds_dU_right = zeros(T, 2)
    t1_right = get_pressure_double_prime(ts, U_right[1]) / (2 * cR)
    ds_dU_right[1] = wR * (t1_right - vR * t2_right)
    ds_dU_right[2] = wR * t2_right
    

    return ds_dU_left, ds_dU_right
end


function numerical_flux(ts::TransientSimulator, U_left::Vector{T}, U_right::Vector{T}, pressure_coeff::Float64)::Vector{T} where {T<:Real}
    
    s = rusanov_speed(ts, U_left, U_right)

    f1 = flux(ts, U_right, pressure_coeff)
    f2 = flux(ts, U_left, pressure_coeff)
    
    return  (f1 + f2 - s * (U_right - U_left) ) / 2.0
end


function jacobian_numerical_flux(ts::TransientSimulator, U_left::Vector{T}, U_right::Vector{T}, pressure_coeff::Float64)::Tuple{ Matrix{T}, Matrix{T} } where {T<:Real}
    J_1 = jacobian_flux(ts, U_left, pressure_coeff)
    J_2 = jacobian_flux(ts, U_right, pressure_coeff)
    s = rusanov_speed(ts, U_left, U_right)
    ds_dU_left, ds_dU_right = rusanov_speed_derivatives(ts, U_left, U_right)

    J_left = (J_1 + s * I - (U_right - U_left) * ds_dU_left') / 2.0
    J_right = (J_2 - s * I - (U_right - U_left) * ds_dU_right') / 2.0
    return J_left, J_right
end

function flux(ts::TransientSimulator, U::Vector{T}, pressure_coeff::Float64)::Vector{T} where {T<:Real}
    inertial_switch = params(ts, :inertial_flag) ? 1 : 0
    var1  = U[2]
    var2  = inertial_switch * (U[2]^2 / U[1]) + pressure_coeff * get_pressure(ts, U[1])
    return [var1, var2]
end

function jacobian_flux(ts::TransientSimulator, U::Vector{T}, pressure_coeff::Float64)::Matrix{T} where {T<:Real}
    inertial_switch = params(ts, :inertial_flag) ? 1 : 0
    J = zeros(T, length(U), length(U))
    J[1, 1] = 0.0
    J[1, 2] = 1.0
    J[2, 1] = pressure_coeff * get_pressure_prime(ts, U[1]) - inertial_switch * (U[2]^2 / (U[1]^2))
    J[2, 2] = inertial_switch * (2.0 * U[2] / U[1])
    return J
end

function source(U::Vector{T}, drag_coeff::Float64, grav_coeff::Float64)::Vector{T} where {T<:Real}
    s = zeros(T, 2)
    s[2] += -drag_coeff * U[2] * abs(U[2]) / U[1]  # friction term
    s[2] += grav_coeff * U[1]  # gravity term 
    return s
end

function jacobian_source(U::Vector{T}, drag_coeff::Float64, grav_coeff::Float64)::Matrix{T} where {T<:Real}
    # Placeholder for Jacobian of source term. This should compute the Jacobian matrix of the source term with respect to U.
    J = zeros(T, length(U), length(U))
    J[2, 1] += grav_coeff # d(gravity term)/d(rho)
    J[2, 1] += drag_coeff * (U[2] * abs(U[2]) / (U[1]^2))  # d(friction term)/d(rho)
    J[2, 2] += -drag_coeff * (2.0 * abs(U[2]) / U[1])  # d(friction term)/d(phi)
    return J
end

# function assemble_local_residual!(ts::TransientSimulator, r::Vector{T}, r_local::Vector{T}, index::Int64, offset::Int64) where {T<:Real}
#     r[index] += r_local[1]
#     r[index + offset] += r_local[2]
#     return
# end

"""
Build Jacobian J = d r / d x for x = [rho; phi], where
r[1:n]     = mass-balance residual
r[n+1:2n]  = momentum-balance residual
matching solve_pipe_state! residual definitions.
"""
function _pipe_jacobian_hyperbolic!(
    J::AbstractMatrix{T},
    ts::TransientSimulator,
    pipe_id::Int64,
    x::Vector{T},
    rho_from::T,
    rho_to::T) where {T<:Real}

    n = div(length(x), 2)
    rho = x[1:n]
    phi = x[n+1:2*n]
    @assert length(phi) == n
    @assert size(J, 1) == 2n && size(J, 2) == 2n

    if J isa SparseMatrixCSC
        fill!(J.nzval, zero(T))
    else
        fill!(J, zero(T))
    end
    I = zeros(T, 2, 2)
    I[1, 1] = T(1.0)
    I[2, 2] = T(1.0)

    dx = ref(ts, :pipe, pipe_id)["dx"]
    dt = params(ts, :dt)
    mu =  dt / dx
    drag_coeff = ref(ts, :pipe, pipe_id, "friction_factor") /
        (2 * ref(ts, :pipe, pipe_id, "diameter"))
    g_coeff = ref(ts, :pipe, pipe_id, "sin_incline")  / (nominal_values(ts, :froude_num))^2
    p_coeff = nominal_values(ts, :euler_num) / (nominal_values(ts, :mach_num))^2

    # Extrapolate edge fluxes.
    phi_from, phi_to = get_bdry_flux_2nd_order_char_extrapolation(ts, rho_from, rho_to, rho[1], rho[2], rho[n-1], rho[n], phi[1], phi[2], phi[n-1], phi[n])
    Bdry_left_mat_1, Bdry_left_mat_2, Bdry_right_mat_n_minus_1, Bdry_right_mat_n = get_bdry_flux_matrices_for_2nd_order_char_extrapolation(ts, rho_from, rho_to, rho[1], rho[2], rho[n-1], rho[n], phi[1], phi[2], phi[n-1], phi[n])
    

    for i = 1 : n
        
        d_mat = zeros(T, 2, 2)
        u_mat = zeros(T, 2, 2)
        l_mat = zeros(T, 2, 2)
        U = Vector{T}([ rho[i], phi[i] ])
        d_mat += I - dt * jacobian_source(U, drag_coeff, g_coeff) 
        if i == 1
            U_right = [ rho[i+1], phi[i+1] ]
            U_left = [rho_from, phi_from]
            J_i_plus_half_left, J_i_plus_half_right = jacobian_numerical_flux(ts, U, U_right, p_coeff)
            J_i_minus_half_left, J_i_minus_half_right = jacobian_numerical_flux(ts, U_left, U, p_coeff)
            d_mat += mu * (J_i_plus_half_left - J_i_minus_half_right -  J_i_minus_half_left * Bdry_left_mat_1)
            u_mat +=  mu * (J_i_plus_half_right -  J_i_minus_half_left * Bdry_left_mat_2)
        elseif i == n
            U_left = [rho[i-1], phi[i-1] ]
            U_right = [rho_to, phi_to]
            J_i_plus_half_left, J_i_plus_half_right = jacobian_numerical_flux(ts, U, U_right, p_coeff)
            J_i_minus_half_left, J_i_minus_half_right = jacobian_numerical_flux(ts, U_left, U, p_coeff)
            d_mat += mu * (J_i_plus_half_left - J_i_minus_half_right +  J_i_plus_half_right * Bdry_right_mat_n)
            l_mat +=  mu * (J_i_plus_half_right  * Bdry_right_mat_n_minus_1 -  J_i_minus_half_left)
        else
            U_left = [ rho[i-1], phi[i-1] ]
            U_right =[ rho[i+1], phi[i+1] ]
            J_i_plus_half_left, J_i_plus_half_right = jacobian_numerical_flux(ts, U, U_right, p_coeff)
            J_i_minus_half_left, J_i_minus_half_right = jacobian_numerical_flux(ts, U_left, U, p_coeff)
            d_mat += mu * (J_i_plus_half_left - J_i_minus_half_right)
            u_mat +=  mu * (J_i_plus_half_right)
            l_mat +=  mu * (-J_i_minus_half_left)     
        end
        assemble_local_jacobian!(ts, J, d_mat, i, i, n)
        if i < n
            assemble_local_jacobian!(ts, J, u_mat, i, i+1, n)
        end
        if i > 1
            assemble_local_jacobian!(ts, J, l_mat, i, i-1, n)
        end
    end
    return J
end

# function assemble_local_jacobian!(ts::TransientSimulator, J::AbstractMatrix{T}, local_mat::AbstractMatrix{T}, row_num::Int64, col_num::Int64, offset::Int64) where {T<:Real}
#     J[row_num, col_num] += local_mat[1, 1]
#     J[row_num, col_num + offset] += local_mat[1, 2]
#     J[row_num + offset, col_num] += local_mat[2, 1]
#     J[row_num + offset, col_num + offset] += local_mat[2, 2]
#     return
# end
