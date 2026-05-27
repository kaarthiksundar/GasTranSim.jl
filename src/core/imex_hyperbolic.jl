

# how to choose timestep so that pipegrid does not change every step and each pipe advances the same  time step
# for each time step, precompute convective flux i+1/2 and i-1/2 ve tors, and pass as argument to pipe residual func
function imex_hyperbolic_step!(ts::TransientSimulator, run_type::Symbol)
    advance_junction_pressures!(ts, :imex_hyperbolic, run_type)
    return
end


function initialize_pipe_grid!(ts::TransientSimulator, ::Val{:imex_hyperbolic})

    dt_factor = 1.0
    while true
        try
            ts.params[:base_dt] = dt_factor * ts.params[:base_dt]
            for (key, pipe) in ref(ts, :pipe)
                # CFL condition (2v)*dt/dx <= 0.9 =>  dx >= (2v)* dt/0.9
                # with nondim dt, dx, we have nondim_dx >= = (2 * v_nondim) * nondim_dt/ 0.9.
                # with dx = L / (n+1) so L / (n+1) >= ( 2v * dt)/0.9
                # or n+1 / L <= 0.9 / ( 2v * dt)
                # or n + 1 <= (0.9 * L) / ( 2v * dt)
                # can see that c instead of 2v is a stronger restriction
                # typical values of v/c are 0.05, 2v/c ub 0.2 suppose
                # 2v/ c = 0.2, (2v/v0) * (v0/c) = 0.2
                # or (2v/v0) = 0.2 / M
                
                c_number = params(ts, :courant_number)
                non_dim_wave_speed = 0.2 / nominal_values(ts, :mach_num)
                num_segments = ( pipe["length"] * c_number ) / ( non_dim_wave_speed * params(ts, :base_dt) ) - 1
                if num_segments < 10
                    throw(CFLException(string(key)))
                end
                num_cvs = floor(Int64, num_segments)
                
                ref(ts, :pipe, key)["num_discretization_points"] = num_cvs
                ref(ts, :pipe, key)["dx"] = pipe["length"] / (num_cvs + 1)
                ref(ts, :pipe, key)["density_profile"] = zeros(Float64, num_cvs)
                ref(ts, :pipe, key)["mass_flux_profile"] = zeros(Float64, num_cvs)
            end
            break
        catch err

            if isa(err, CFLException)
                dt_factor *= 0.5
                println("CFL condition failed. Reducing Δt by half")
            else
                rethrow(err)
            end

        end

    
        
    end
    return
end

function initialize_pipe_state!(ts::TransientSimulator, ::Val{:imex_hyperbolic})
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



function _solve_pipe_state_imex_hyperbolic!(
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


    residual_fun! = (r, x) -> _pipe_residual_imex_hyperbolic!(
        r,
        ts,
        pipe_id,
        x,
        rho_old,
        phi_old,
        rho_from,
        rho_to
    )
    Jacobian_fun! = (J, x) -> _pipe_jacobian_imex_hyperbolic!(
        J,
        ts,
        pipe_id,
        x,
        rho_from,
        rho_to
    )


    x, converged, _, res_norm = solve_newton_basic!(x, residual_fun!, Jacobian_fun!, tol = 1e-6, max_iter = 100)

    
    converged || throw(DomainError(res_norm, "Newton solver did not converge for pipe $pipe_id"))

    # phi_from, phi_to = get_bry_flux_imex_char_extrapolation(ts, rho_from, rho_to, x[1], x[n], x[n+1], x[n+n])
    phi_from, phi_to = get_bdry_flux_imex_2nd_order_char_extrapolation(ts, rho_from, rho_to, x[1], x[2], x[n-1], x[n], x[n+1], x[n+2], x[n+n-1], x[n+n])

    if T == Float64
        ref(ts, :pipe, pipe_id)["rho"] = x[1:n]
        ref(ts, :pipe, pipe_id)["phi"] = x[n+1:2*n]
        ref(ts, :pipe, pipe_id)["fr_mass_flux"] = phi_from
        ref(ts, :pipe, pipe_id)["to_mass_flux"] = phi_to
    end


    end_flows = T[area * phi_from, area * phi_to]
    
    return end_flows

end

function imex_characteristic_potential(ts::TransientSimulator, rho::T)::T where {T<:Real}
    nquad = 100
    delta_rho = rho / T(nquad)
    acc = zero(T)
    # use mid pt rule to avoid  zero density
    for i in 0:(nquad - 1)
        rho_i = (T(i) + T(0.5)) * delta_rho
        integrand = sqrt(get_pressure_prime(ts, rho_i))
        acc += integrand
    end
    return acc * delta_rho
end

function imex_characteristic_potential_prime(ts::TransientSimulator, rho::T)::T where {T<:Real}
    derivative =  sqrt(get_pressure_prime(ts, rho))
    return derivative
end


function get_bry_flux_imex_char_extrapolation(ts::TransientSimulator, rho_from::T, rho_to::T, rho_1::T, rho_n::T, phi_1::T, phi_n::T)::Tuple{T, T} where {T<:Real}

    
    pot_from = imex_characteristic_potential(ts, rho_from)
    pot_1 = imex_characteristic_potential(ts, rho_1)
    pot_to = imex_characteristic_potential(ts, rho_to)
    pot_n = imex_characteristic_potential(ts, rho_n)

    phi_from = phi_1  - pot_1 + pot_from
    phi_to = phi_n   + pot_n - pot_to

    return phi_from, phi_to
end

function get_bry_flux_imex_matrices_for_char_extrapolation(ts::TransientSimulator, rho_from::T, rho_to::T, rho_1::T, rho_n::T,phi_1::T, phi_n::T)::Tuple{Matrix{T}, Matrix{T}} where {T<:Real}

    
    Bdry_left_mat = [T(0) T(0); -imex_characteristic_potential_prime(ts, rho_1)  T(1)]
    Bdry_right_mat = [T(0) T(0);  imex_characteristic_potential_prime(ts, rho_n)  T(1)]
    
    return Bdry_left_mat, Bdry_right_mat
end

function get_bdry_flux_imex_2nd_order_char_extrapolation(ts::TransientSimulator, rho_from::T, rho_to::T, rho_1::T, rho_2::T, rho_n_minus_1::T, rho_n::T, phi_1::T, phi_2::T, phi_n_minus_1::T, phi_n::T)::Tuple{T, T} where {T<:Real}

    pot_from = imex_characteristic_potential(ts, rho_from)
    pot_1 = imex_characteristic_potential(ts, rho_1)
    pot_2 = imex_characteristic_potential(ts, rho_2)

    pot_to = imex_characteristic_potential(ts, rho_to)
    pot_n = imex_characteristic_potential(ts, rho_n)
    pot_n_minus_1 = imex_characteristic_potential(ts, rho_n_minus_1)

    

    phi_from = (3/2) * (phi_1   - pot_1) - (1/2) * (phi_2  - pot_2) + pot_from 
    phi_to = (3/2) * (phi_n  + pot_n) - (1/2) * (phi_n_minus_1  + pot_n_minus_1) - pot_to

    

    return phi_from, phi_to
end

function get_bdry_flux_imex_matrices_for_2nd_order_char_extrapolation(ts::TransientSimulator, rho_from::T, rho_to::T, rho_1::T, rho_2::T, rho_n_minus_1::T, rho_n::T, phi_1::T, phi_2::T, phi_n_minus_1::T, phi_n::T)::Tuple{Matrix{T}, Matrix{T}, Matrix{T}, Matrix{T}} where {T<:Real}

    

    Bdry_left_mat_1 = [
        T(0) T(0);
        -T(3 / 2) * imex_characteristic_potential_prime(ts, rho_1) T(3 / 2)
    ]
    Bdry_left_mat_2 = [
        T(0) T(0);
        T(1 / 2) * imex_characteristic_potential_prime(ts, rho_2) -T(1 / 2)
    ]
    Bdry_right_mat_n_minus_1 = [
        T(0) T(0);
        -T(1 / 2) * imex_characteristic_potential_prime(ts, rho_n_minus_1) -T(1 / 2)
    ]
    Bdry_right_mat_n = [
        T(0) T(0);
        T(3 / 2) * imex_characteristic_potential_prime(ts, rho_n) T(3 / 2)
    ]
    

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
function _pipe_residual_imex_hyperbolic!(
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

    

    phi_from, phi_to = get_bdry_flux_imex_2nd_order_char_extrapolation(ts, rho_from, rho_to, rho[1], rho[2], rho[n-1], rho[n], phi[1], phi[2], phi[n-1], phi[n])


    for i = 1 : n
        fill!(r_local, zero(T))

        U =  [ rho[i], phi[i] ]
        U_old = [rho_old[i], phi_old[i]]
        r_local += U - U_old 
        r_local += - dt * imex_source(U, drag_coeff, g_coeff)
        if i == 1
            flux_i_plus_half  = imex_numerical_flux(ts, U, [rho[i+1], phi[i+1]], p_coeff) 
            conv_i_plus_half = numerical_convection_flux(ts, pipe_id, rho_old, phi_old, 1)
            flux_i_minus_half = imex_numerical_flux(ts, [rho_from, phi_from], U, p_coeff) 
            conv_i_minus_half = numerical_convection_flux(ts, pipe_id,rho_old, phi_old, 0)
        elseif i == n
            flux_i_plus_half  = imex_numerical_flux(ts, U, [rho_to, phi_to], p_coeff) 
            conv_i_plus_half = numerical_convection_flux(ts, pipe_id,rho_old, phi_old, n)
            flux_i_minus_half = imex_numerical_flux(ts, [rho[i-1], phi[i-1]], U, p_coeff)
            conv_i_minus_half = numerical_convection_flux(ts, pipe_id,rho_old, phi_old, n-1)    
        else
            flux_i_plus_half  = imex_numerical_flux(ts, U, [ rho[i+1], phi[i+1] ], p_coeff)
            conv_i_plus_half = numerical_convection_flux(ts, pipe_id,rho_old, phi_old, i) 
            flux_i_minus_half = imex_numerical_flux(ts, [rho[i-1], phi[i-1]], U, p_coeff)
            conv_i_minus_half = numerical_convection_flux(ts, pipe_id,rho_old, phi_old, i-1)
        end
        r_local +=  mu * ( flux_i_plus_half  - flux_i_minus_half  )  + mu * (conv_i_plus_half   - conv_i_minus_half )
        assemble_local_residual!(ts, r, r_local, i, n)
    end
    
    return r
end


function imex_rusanov_speed(ts::TransientSimulator, U_left::AbstractVector{T}, U_right::AbstractVector{T})::T where {T<:Real}
    epsilon = T(1e-6)
    cL = sqrt(get_pressure_prime(ts, U_left[1])) # speed of sound on the left
    cR = sqrt(get_pressure_prime(ts, U_right[1])) # speed of sound on the right
    return smooth_max(cL, cR, epsilon)
end



function imex_rusanov_speed_derivatives(ts::TransientSimulator, U_left::AbstractVector{T}, U_right::AbstractVector{T})::Tuple{Vector{T}, Vector{T}} where {T<:Real}
    
    epsilon = T(1e-6)
    cL = sqrt(get_pressure_prime(ts, U_left[1])) # speed of sound on the left
    cR = sqrt(get_pressure_prime(ts, U_right[1])) # speed of sound on the right

    sL = cL
    sR = cR
    
    wL  = T(1/2) * (1 + (sL - sR) / sqrt((sL - sR)^2 + epsilon^2) )
    wR  = T(1/2) * (1 - (sL - sR) / sqrt((sL - sR)^2 + epsilon^2) )

    ds_dU_left = zeros(T, 2)

    t1_left = get_pressure_double_prime(ts, U_left[1]) / (2 * cL)
    ds_dU_left[1] = wL * (t1_left)
    ds_dU_left[2] = T(0)

    ds_dU_right = zeros(T, 2)
    t1_right = get_pressure_double_prime(ts, U_right[1]) / (2 * cR)
    ds_dU_right[1] = wR * (t1_right)
    ds_dU_right[2] = T(0)
    

    return ds_dU_left, ds_dU_right
end


function imex_numerical_flux(ts::TransientSimulator, U_left::Vector{T}, U_right::Vector{T}, pressure_coeff::Float64)::Vector{T} where {T<:Real}
    
    s = imex_rusanov_speed(ts, U_left, U_right)

    f1 = imex_flux(ts, U_right, pressure_coeff)
    f2 = imex_flux(ts, U_left, pressure_coeff)
    
    return  (f1 + f2 - s * (U_right - U_left) ) / 2.0
end


function imex_jacobian_numerical_flux(ts::TransientSimulator, U_left::Vector{T}, U_right::Vector{T}, pressure_coeff::Float64)::Tuple{ Matrix{T}, Matrix{T} } where {T<:Real}
    J_1 = imex_jacobian_flux(ts, U_left, pressure_coeff)
    J_2 = imex_jacobian_flux(ts, U_right, pressure_coeff)
    s = imex_rusanov_speed(ts, U_left, U_right)
    ds_dU_left, ds_dU_right = imex_rusanov_speed_derivatives(ts, U_left, U_right)

    J_left = (J_1 + s * I - (U_right - U_left) * ds_dU_left') / 2.0
    J_right = (J_2 - s * I - (U_right - U_left) * ds_dU_right') / 2.0
    return J_left, J_right
end

function imex_flux(ts::TransientSimulator, U::Vector{T}, pressure_coeff::Float64)::Vector{T} where {T<:Real}
    var1  = U[2]
    var2  = pressure_coeff * get_pressure(ts, U[1])
    return [var1, var2]
end

function imex_jacobian_flux(ts::TransientSimulator, U::Vector{T}, pressure_coeff::Float64)::Matrix{T} where {T<:Real}
    J = zeros(T, length(U), length(U))
    J[1, 1] = T(0)
    J[1, 2] = T(1)
    J[2, 1] = pressure_coeff * get_pressure_prime(ts, U[1]) 
    return J
end

function convection_flux_explicit(ts::TransientSimulator, U::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    return [T(0), U[2]^2 / U[1]]
end

function convection_rusanov_speed(ts::TransientSimulator, U_left::AbstractVector{T}, U_right::AbstractVector{T})::T where {T<:Real}
    epsilon = T(1e-6)
    vL = U_left[2] / U_left[1] # velocity on the left
    vR = U_right[2] / U_right[1] # velocity on the right
    sL = smooth_abs(2 * vL, epsilon)
    sR = smooth_abs(2 * vR, epsilon)
    return smooth_max(sL, sR, epsilon)
end

function numerical_convection_flux(ts::TransientSimulator, pipe_id::Integer, rho_old::AbstractVector{T},
    phi_old::AbstractVector{T}, i_left::Integer)::AbstractVector{T} where {T<:Real}
    
    if i_left == 0
        fr_node = ref(ts, :pipe, pipe_id, "fr_node")
        rho_fr = T(get_density(ts, ref(ts, :node, fr_node)["pressure"]))
        phi_fr = T(ref(ts, :pipe, pipe_id)["fr_mass_flux"])
        U_left = [rho_fr, phi_fr]
        U_right = [ rho_old[i_left + 1], phi_old[i_left + 1] ]

    elseif i_left == length(rho_old)
        to_node = ref(ts, :pipe, pipe_id, "to_node")
        rho_to = T(get_density(ts, ref(ts, :node, to_node)["pressure"]))
        phi_to = T(ref(ts, :pipe, pipe_id)["to_mass_flux"])
        U_right = [rho_to, phi_to]
        U_left = [ rho_old[i_left], phi_old[i_left] ]
    else
        U_right = [ rho_old[i_left + 1], phi_old[i_left + 1] ]
        U_left = [ rho_old[i_left], phi_old[i_left] ]
    end
    
    s = convection_rusanov_speed(ts, U_left, U_right)

    f1 = convection_flux_explicit(ts, U_right)
    f2 = convection_flux_explicit(ts, U_left)
    
    return  (f1 + f2 - s * (U_right - U_left) ) / 2.0
end



function imex_source(U::Vector{T}, drag_coeff::Float64, grav_coeff::Float64)::Vector{T} where {T<:Real}
    s = zeros(T, 2)
    s[2] += -drag_coeff * U[2] * abs(U[2]) / U[1]  # friction term
    s[2] += grav_coeff * U[1]  # gravity term 
    return s
end

function imex_jacobian_source(U::Vector{T}, drag_coeff::Float64, grav_coeff::Float64)::Matrix{T} where {T<:Real}
    # Placeholder for Jacobian of source term. This should compute the Jacobian matrix of the source term with respect to U.
    J = zeros(T, length(U), length(U))
    J[2, 1] += grav_coeff # d(gravity term)/d(rho)
    J[2, 1] += drag_coeff * (U[2] * abs(U[2]) / (U[1]^2))  # d(friction term)/d(rho)
    J[2, 2] += -drag_coeff * (2.0 * abs(U[2]) / U[1])  # d(friction term)/d(phi)
    return J
end



"""
Build Jacobian J = d r / d x for x = [rho; phi], where
r[1:n]     = mass-balance residual
r[n+1:2n]  = momentum-balance residual
matching solve_pipe_state! residual definitions.
"""
function _pipe_jacobian_imex_hyperbolic!(
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
    phi_from, phi_to = get_bdry_flux_imex_2nd_order_char_extrapolation(ts, rho_from, rho_to, rho[1], rho[2], rho[n-1], rho[n], phi[1], phi[2], phi[n-1], phi[n])
    Bdry_left_mat_1, Bdry_left_mat_2, Bdry_right_mat_n_minus_1, Bdry_right_mat_n = get_bdry_flux_imex_matrices_for_2nd_order_char_extrapolation(ts, rho_from, rho_to, rho[1], rho[2], rho[n-1], rho[n], phi[1], phi[2], phi[n-1], phi[n])
    

    for i = 1 : n
        
        d_mat = zeros(T, 2, 2)
        u_mat = zeros(T, 2, 2)
        l_mat = zeros(T, 2, 2)
        U = Vector{T}([ rho[i], phi[i] ])
        d_mat += I - dt * imex_jacobian_source(U, drag_coeff, g_coeff) 
        if i == 1
            U_right = [ rho[i+1], phi[i+1] ]
            U_left = [rho_from, phi_from]
            J_i_plus_half_left, J_i_plus_half_right = imex_jacobian_numerical_flux(ts, U, U_right, p_coeff)
            J_i_minus_half_left, J_i_minus_half_right = imex_jacobian_numerical_flux(ts, U_left, U, p_coeff)
            d_mat += mu * (J_i_plus_half_left - J_i_minus_half_right -  J_i_minus_half_left * Bdry_left_mat_1)
            u_mat +=  mu * (J_i_plus_half_right -  J_i_minus_half_left * Bdry_left_mat_2)
        elseif i == n
            U_left = [rho[i-1], phi[i-1] ]
            U_right = [rho_to, phi_to]
            J_i_plus_half_left, J_i_plus_half_right = imex_jacobian_numerical_flux(ts, U, U_right, p_coeff)
            J_i_minus_half_left, J_i_minus_half_right = imex_jacobian_numerical_flux(ts, U_left, U, p_coeff)
            d_mat += mu * (J_i_plus_half_left - J_i_minus_half_right +  J_i_plus_half_right * Bdry_right_mat_n)
            l_mat +=  mu * (J_i_plus_half_right  * Bdry_right_mat_n_minus_1 -  J_i_minus_half_left)
        else
            U_left = [ rho[i-1], phi[i-1] ]
            U_right =[ rho[i+1], phi[i+1] ]
            J_i_plus_half_left, J_i_plus_half_right = imex_jacobian_numerical_flux(ts, U, U_right, p_coeff)
            J_i_minus_half_left, J_i_minus_half_right = imex_jacobian_numerical_flux(ts, U_left, U, p_coeff)
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


