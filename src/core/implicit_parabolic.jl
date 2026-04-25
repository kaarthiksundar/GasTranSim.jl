function implicit_parabolic_step!(ts::TransientSimulator, run_type::Symbol)
    implicit_advance_junction_pressures!(ts, run_type) # pressure (n+1), flux (n+1/2)
    # update density and mass_flux profile
    # update end fluxes
    return
end


function initialize_pipe_grid!(ts::TransientSimulator, ::Val{:implicit_parabolic})
    for (key, pipe) in ref(ts, :pipe)
        # CFL condition c*dt/dx <= 0.9 => dx >= c*dt/0.9
        # with nondim dt, dx, we have nondim_dt/ nondim_dx < = 0.9 * mach_no
        
        # c_inv = nominal_values(ts, :mach_num)
        # num_segments =
        #     c_inv * (pipe["length"] * params(ts, :courant_number)) / params(ts, :base_dt)
        # if num_segments < 1
        #     throw(CFLException(string(key)))
        # end
        # n = floor(Int64, num_segments) + 1
        n= 61
        ref(ts, :pipe, key)["num_discretization_points"] = n
        ref(ts, :pipe, key)["dx"] = pipe["length"] / (n - 1)
        ref(ts, :pipe, key)["density_profile"] = zeros(Float64, n)
        ref(ts, :pipe, key)["mass_flux_profile"] = zeros(Float64, n)
    end
    return
end

function initialize_pipe_state!(ts::TransientSimulator, ::Val{:implicit_parabolic})
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
        if is_steady
            initial_mass_flux = initial_pipe_mass_flow(ts, key)(0.0) / area
            fill!(pipe["mass_flux_profile"], initial_mass_flux)
            initial_fr_pressure = ref(ts, :node, fr_node, "pressure")
            initial_to_pressure = ref(ts, :node, to_node, "pressure")
            density_at_first_sq = get_density(ts, initial_fr_pressure)^2
            density_at_last_sq = get_density(ts, initial_to_pressure)^2
            dL = dx / L
            pipe["density_profile"][1:n] = [
                sqrt(
                    density_at_last_sq * (i - 1) * dL + density_at_first_sq * (n - i) * dL,
                ) for i = 1:n
            ]
            pipe["fr_minus_mass_flux"] = initial_mass_flux # dx
            pipe["to_minus_mass_flux"] = initial_mass_flux # L-dx
            pipe["fr_mass_flux"] = initial_mass_flux # 0
            pipe["to_mass_flux"] = initial_mass_flux # L
        else
            flow_spl = initial_pipe_mass_flow(ts, key)
            pressure_spl = initial_pipe_pressure(ts, key)
            x_arr = LinRange(0, L, n)
            get_coeffs(flow_spl)[1]
            pipe["mass_flux_profile"][1:n] =
                [flow_spl(x) for x in x_arr] ./ area
            pipe["fr_minus_mass_flux"] = pipe["mass_flux_profile"][2]
            pipe["to_minus_mass_flux"] = pipe["mass_flux_profile"][n-1]
            pipe["fr_mass_flux"] = pipe["mass_flux_profile"][1]
            pipe["to_mass_flux"] = pipe["mass_flux_profile"][end]
            pipe["density_profile"][1:n] = [get_density(ts, pressure_spl(x)) for x in x_arr]
        end
        pipe["phi"] = zeros(n)
        pipe["rho"] = zeros(n)
    end
    return
end

function solve_pipe_state!(
    ts::TransientSimulator,
    pipe_id::Int64,
    rho_from::T,
    rho_to::T,
)::Vector{T} where {T<:Real}
    n = ref(ts, :pipe, pipe_id)["num_discretization_points"]
    x = zeros(T, 2 * n)
    x[1:n] = T.(ref(ts, :pipe, pipe_id)["density_profile"])
    x[n+1:2*n] = T.(ref(ts, :pipe, pipe_id)["mass_flux_profile"])
    rho_old = copy(x[1:n])
    area = T(ref(ts, :pipe, pipe_id)["area"])

    residual_fun! = (r, x) -> _pipe_residual!(
        r,
        ts,
        pipe_id,
        x,
        rho_old,
        rho_from,
        rho_to;
        inertial_flag = zero(T),
    )
    Jacobian_fun! = (J, x) -> _pipe_jacobian!(
        J,
        ts,
        pipe_id,
        x;
        inertial_flag = zero(T),
    )


    x, converged, _, res_norm = solve_newton_basic!(x, residual_fun!, Jacobian_fun!)

    
    converged || throw(DomainError(res_norm, "Newton solver did not converge for pipe $pipe_id"))

    if T == Float64
        ref(ts, :pipe, pipe_id)["rho"] = x[1:n]
        ref(ts, :pipe, pipe_id)["phi"] = x[n+1:2*n]
    end

    end_flows = T[area * x[n+1], area * x[2*n]]
    
    return end_flows

end

"""
Residual for implicit parabolic pipe solve, matching solve_pipe_state! layout.

Unknown/state ordering:
- x[1:n]     -> rho
- x[n+1:2n]  -> phi

Residual ordering:
- r[1:n]     -> mass equations
- r[n+1:2n]  -> momentum equations
"""
function _pipe_residual!(
    r::Vector{T},
    ts::TransientSimulator,
    pipe_id::Int64,
    x::Vector{T},
    rho_old::Vector{T},
    rho_from::T,
    rho_to::T;
    inertial_flag::T) where {T<:Real}

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

    # -------------------------
    # Mass balance residual
    # -------------------------
    r[1] = rho[1] - rho_from
    r[2:n-1] .= rho[2:n-1] .- rho_old[2:n-1] .+ (dt / dx) .* (phi[2:n-1] .- phi[1:n-2])
    r[n] = rho[n] - rho_to

    # -------------------------
    # Momentum balance residual
    # -------------------------
    rho_x = zeros(T, n)
    rho_x[1] = (rho[2] - rho[1]) / dx
    rho_x[2:n-1] .= (rho[3:n] .- rho[1:n-2]) ./ (2.0 * dx)
    rho_x[n] = (rho[n] - rho[n-1]) / dx

    phi_sq = phi .^ 2
    conv_term = zeros(T, n)
    conv_term[1] = (phi_sq[2] / rho[2] - phi_sq[1] / rho[1]) / dx
    conv_term[2:n-1] .= (phi_sq[3:n] ./ rho[3:n] .- phi_sq[1:n-2] ./ rho[1:n-2]) ./ (2.0 * dx)
    conv_term[n] = (phi_sq[n] / rho[n] - phi_sq[n-1] / rho[n-1]) / dx

    c2 = get_pressure_prime.(Ref(ts), rho)  # vectorized in your types.jl
    beta = ref(ts, :pipe, pipe_id, "friction_factor") /
        (2 * ref(ts, :pipe, pipe_id, "diameter"))
    nondim = nominal_values(ts, :euler_num) / (nominal_values(ts, :mach_num))^2
    phi_fun = phi .* abs.(phi)
    

    r[n+1:2n] .=
        nondim .* c2 .* rho_x .+
        inertial_flag .* conv_term .+
        beta .* phi_fun ./ rho .-
        ref(ts, :pipe, pipe_id, "sin_incline") .* rho / (nominal_values(ts, :froude_num))^2

    return r
end

# Convenience allocating version
function _pipe_residual(
    ts::TransientSimulator,
    pipe_id::Int64,
    rho::Vector{Float64},
    phi::Vector{Float64},
    rho_old::Vector{Float64},
    rho_from::Float64,
    rho_to::Float64;
    kwargs...)
    r = zeros(Float64, 2length(rho))
    return pipe_residual!(r, ts, pipe_id, rho, phi, rho_old, rho_from, rho_to; kwargs...)
end


"""
Build Jacobian J = d r / d x for x = [rho; phi], where
r[1:n]     = mass-balance residual
r[n+1:2n]  = momentum-balance residual
matching solve_pipe_state! residual definitions.
"""
function _pipe_jacobian!(
    J::AbstractMatrix{T},
    ts::TransientSimulator,
    pipe_id::Int64,
    x::Vector{T};
    inertial_flag::T) where {T<:Real}

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

    dx = ref(ts, :pipe, pipe_id)["dx"]

    rho_col(i) = i
    phi_col(i) = n + i
    mom_row(i) = n + i

    # -------------------------
    # Mass rows: r[1:n]
    # -------------------------
    J[1, rho_col(1)] = 1.0
    for i in 2:n-1
        J[i, rho_col(i)] = 1.0
        J[i, phi_col(i)] = 1.0 / dx * params(ts, :dt)
        J[i, phi_col(i - 1)] = -1.0 / dx * params(ts, :dt)
    end
    J[n, rho_col(n)] = 1.0

    # c2 = dp/drho, c2p = d2p/drho2
    c2 = get_pressure_prime.(Ref(ts), rho)
    c2p = get_pressure_double_prime.(Ref(ts), rho)
    nondim = nominal_values(ts, :euler_num) / (nominal_values(ts, :mach_num))^2
    


    # helper for conv-term contributions from g_k = phi_k^2 / rho_k
    function add_conv_g!(row::Int, k::Int, coeff::Real)
        J[row, rho_col(k)] += inertial_flag * coeff * (-(phi[k]^2) / (rho[k]^2))
        J[row, phi_col(k)] += inertial_flag * coeff * (2.0 * phi[k] / rho[k])
    end

    for i in 1:n
        row = mom_row(i)

        # rho_x(i) and its rho-derivatives
        rho_x_i = zero(T)
        if i == 1
            rho_x_i = (rho[2] - rho[1]) / dx
            J[row, rho_col(1)] += nondim * c2[i] * (-1.0 / dx)
            J[row, rho_col(2)] += nondim * c2[i] * ( 1.0 / dx)
        elseif i == n
            rho_x_i = (rho[n] - rho[n - 1]) / dx
            J[row, rho_col(n - 1)] += nondim * c2[i] * (-1.0 / dx)
            J[row, rho_col(n)]     += nondim * c2[i] * ( 1.0 / dx)
        else
            rho_x_i = (rho[i + 1] - rho[i - 1]) / (2.0 * dx)
            J[row, rho_col(i - 1)] += nondim * c2[i] * (-1.0 / (2.0 * dx))
            J[row, rho_col(i + 1)] += nondim * c2[i] * ( 1.0 / (2.0 * dx))
        end

        # chain term from c2(rho_i) * rho_x(i)
        J[row, rho_col(i)] += nondim * c2p[i] * rho_x_i

        # inertial conv_term
        if i == 1
            add_conv_g!(row, 2,  1.0 / dx)
            add_conv_g!(row, 1, -1.0 / dx)
        elseif i == n
            add_conv_g!(row, n,      1.0 / dx)
            add_conv_g!(row, n - 1, -1.0 / dx)
        else
            add_conv_g!(row, i + 1,  1.0 / (2.0 * dx))
            add_conv_g!(row, i - 1, -1.0 / (2.0 * dx))
        end

        # friction term: beta * phi_i * abs(phi_i) / rho_i
        beta = ref(ts, :pipe, pipe_id, "friction_factor") /
        (2 * ref(ts, :pipe, pipe_id, "diameter"))
        J[row, rho_col(i)] += -beta * phi[i] * abs(phi[i]) / (rho[i]^2)
        J[row, phi_col(i)] +=  beta * (2.0 * abs(phi[i])) / rho[i]  # subgradient at 0 -> 0

        # gravity term: -G * rho_i
        J[row, rho_col(i)] += -ref(ts, :pipe, pipe_id, "sin_incline") / (nominal_values(ts, :froude_num))^2
    end

    return J
end


function implicit_advance_junction_pressures!(ts::TransientSimulator, _run_type::Symbol)
    x_node = get_density.(Ref(ts), form_nodal_pressure_vector(ts))
    
    problem_fun! = (r, J, x) -> assemble_network_problem!(ts, x, r, J)

    x_node, converged, iter, res_norm = NR_solve!(x_node, problem_fun!)
    # println(iter)
    converged || throw(DomainError(res_norm, "Newton solver did not converge for nodal densities"))

    check_limits(ts, x_node)

    # Update nodal pressures in ts.ref using the density solution x_node.
    for node_id = 1:length(x_node)
        p_val = get_pressure(ts, x_node[node_id])
        ref(ts, :node, node_id)["pressure_previous"] = ref(ts, :node, node_id)["pressure"]
        ref(ts, :node, node_id)["pressure"] = p_val
        ref(ts, :node, node_id)["is_updated"] = true
    end

    # update mass flux and density profiles in pipes
    for (pipe_id, pipe) in ref(ts, :pipe)
        pipe["density_profile"] = pipe["rho"]
        pipe["mass_flux_profile"] = pipe["phi"]
        pipe["fr_mass_flux"] = pipe["phi"][1]
        pipe["to_mass_flux"] = pipe["phi"][end]
    end

    return
end

function assemble_network_problem!(ts::TransientSimulator, x::Vector{Float64}, r::Vector{Float64}, J::SparseMatrixCSC{Float64,Int64})

    _assemble_for_nodes_WITHOUT_eqn_nos!(ts, x, r, J)
    _assemble_for_nodes_WITH_eqn_nos!(ts, x, r, J)
    _assemble_all_pipes!(ts, x, r, J)

    return
end



function _assemble_all_pipes!(ts::TransientSimulator, x::Vector{Float64}, r::Vector{Float64}, J::SparseMatrixCSC{Float64,Int64})

    for (pipe_id, pipe) in ref(ts, :pipe)
        to_node = pipe["to_node"]
        fr_node = pipe["fr_node"]
        rho_to = x[to_node]
        rho_from = x[fr_node]

        end_flow_func = bdry_vals -> solve_pipe_state!(ts, pipe_id,bdry_vals[1], bdry_vals[2])
        end_flows = end_flow_func([rho_from, rho_to])
        sensitivity_mat = ForwardDiff.jacobian(end_flow_func, [rho_from, rho_to])

        _assemble_pipe_solve_results!(ts, fr_node, to_node, end_flows, sensitivity_mat, r, J)
    end

    return
end

function _assemble_pipe_solve_results!(ts::TransientSimulator, fr_node::Int64, to_node::Int64, end_flows::Vector{Float64}, sensitivity_mat::AbstractArray,r::Vector{Float64}, J::SparseMatrixCSC{Float64,Int64})
    flow_at_fr_end, flow_at_to_end = end_flows
    s_from_from, s_from_to = sensitivity_mat[1,:]
    s_to_from, s_to_to = sensitivity_mat[2,:]
        # now assemble these into r and J
    eq_num_fr =  ref(ts, :node, fr_node)["eqn_number"]
    if !isnan(eq_num_fr) && ref(ts, :node, fr_node)["is_slack"] == 0
        r[eq_num_fr] += -flow_at_fr_end # positive dir serves as withdrawal from end
        J[eq_num_fr, fr_node] += -s_from_from
        J[eq_num_fr, to_node] += -s_from_to
    end

    eq_num_to =  ref(ts, :node, to_node)["eqn_number"]
    if !isnan(eq_num_to) && ref(ts, :node, to_node)["is_slack"] == 0
        r[eq_num_to] += flow_at_to_end # positive dir serves as injection into end
        J[eq_num_to, fr_node] += s_to_from
        J[eq_num_to, to_node] += s_to_to
    end
    return
end


function _assemble_for_nodes_WITHOUT_eqn_nos!(ts::TransientSimulator, x_node::Vector{Float64}, residual_node::Vector{Float64}, J::SparseMatrixCSC{Float64,Int64})
    for (ci, compressor) in  get(ref(ts), :compressor, Dict())

        to_node = compressor["to_node"]
        from_node = compressor["fr_node"]
        ctrl_type, ctrl_val = control(ts, :compressor, ci, ref(ts, :current_time))
        if ctrl_type == c_ratio_control
            residual_node[ci] = get_pressure(ts, x_node[to_node]) - ctrl_val * get_pressure(ts, x_node[from_node]) 
            J[ci, from_node] = -ctrl_val * get_pressure_prime(ts, x_node[from_node])
            J[ci, to_node]  = get_pressure_prime(ts, x_node[to_node]) 
        elseif ctrl_type == discharge_pressure_control
            residual_node[ci] = x_node[to_node] - get_density(ts, ctrl_val)
            J[ci, to_node] = 1.0
        elseif ctrl_type == flow_control
            # no contribution to Jacobian since flow is known
            continue
        end
    end
    return
end

    

function _assemble_for_nodes_WITH_eqn_nos!(ts::TransientSimulator, x_node::Vector{Float64}, residual_node::Vector{Float64}, J::SparseMatrixCSC{Float64,Int64})
    for (node_id, node) in ref(ts, :node)

        eqn_num = ref(ts, :node, node_id)["eqn_number"]
        if  isnan(eqn_num)
            continue
        end
        ctrl_type, ctrl_val = control(ts, :node, node_id, ref(ts, :current_time))
        if ctrl_type == pressure_control
            residual_node[eqn_num] = x_node[node_id]- get_density(ts, ctrl_val)
            J[eqn_num, node_id] = 1.0
            continue
        elseif ctrl_type == flow_control
            
            rhs_compressor_term = _assemble_for_flow_control_compressors(node_id, ts)
            residual_node[eqn_num] += -ctrl_val + rhs_compressor_term  #ctrl_val is  a withdrawal
            continue
        end
    end
    return
end

function _assemble_for_flow_control_compressors(
    node_id::Int64,
    ts::TransientSimulator)::Real
    out_c = ref(ts, :outgoing_compressors, node_id)
    in_c = ref(ts, :incoming_compressors, node_id)
    rhs_term = 0.0
    t = ref(ts, :current_time)
    for ci in in_c
        ctr, cmpr_val = control(ts, :compressor, ci, t)
        if ctr == flow_control
            rhs_term += cmpr_val # inflow positive
        end
    end
    for co in out_c
        ctr, cmpr_val = control(ts, :compressor, co, t)
        if ctr == flow_control
            rhs_term += (-1.0 * cmpr_val) # outflow negative
        end
    end
    return rhs_term
end



function NR_solve!(
    x::Vector{Float64},
    problem_fun!::Function;
    tol::Float64 = 1e-6,
    max_iter::Int = 100)::Tuple{Vector{Float64},Bool,Int,Float64}

    n = length(x)
    residual = zeros(Float64, n)
    J = spzeros(n, n)

    for iter = 1:max_iter
        fill!(residual, 0.0)
        fill!(J.nzval, 0.0)

        problem_fun!(residual, J, x)
        res_norm = maximum(abs, residual)
        # println(iter, ":", res_norm)
        if res_norm <= tol
            return x, true, iter, res_norm
        end

        delta_x = J \ (-residual)
        x .+= delta_x

        step_norm = maximum(abs, delta_x)
        if step_norm <= tol
            return x, true, iter, res_norm
        end
    end

    fill!(residual, 0.0)
    problem_fun!(residual, J, x)
    return x, false, max_iter, maximum(abs, residual)
end