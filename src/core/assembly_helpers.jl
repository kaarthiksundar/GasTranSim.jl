
function precompute_block!(ts::TransientSimulator, precompute_function::Function)

    if !isnothing(precompute_function)
        precompute_function(ts)
    end

    return
end


function advance_junction_pressures!(ts::TransientSimulator, method::Symbol,  _run_type::Symbol)
    x_node = get_density.(Ref(ts), form_nodal_pressure_vector(ts))
    
    problem_fun! = (r, J, x) -> assemble_network_problem!(ts, method, x, r, J)
    
    x_node, converged, iter, res_norm = NR_solve!(x_node, problem_fun!)
    
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

function assemble_network_problem!(ts::TransientSimulator, method::Symbol, x::Vector{Float64}, r::Vector{Float64}, J::SparseMatrixCSC{Float64,Int64})

    _assemble_for_nodes_WITHOUT_eqn_nos!(ts, x, r, J)
    _assemble_for_nodes_WITH_eqn_nos!(ts, x, r, J)
    _assemble_all_pipes!(ts, method, x, r, J)

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


function _assemble_all_pipes!(ts::TransientSimulator, method::Symbol,  x::Vector{Float64}, r::Vector{Float64}, J::SparseMatrixCSC{Float64,Int64})

    for (pipe_id, pipe) in ref(ts, :pipe)
        to_node = pipe["to_node"]
        fr_node = pipe["fr_node"]
        rho_to = x[to_node]
        rho_from = x[fr_node]
        dx = ref(ts, :pipe, pipe_id)["dx"]
        dt = params(ts, :dt)
        area = ref(ts, :pipe, pipe_id)["area"]

        mu = area * dx / dt

        end_flow_func = bdry_vals -> solve_pipe_state!(ts, method, pipe_id,bdry_vals[1], bdry_vals[2])
        end_flows = end_flow_func([rho_from, rho_to])
        if method == :explicit_hyperbolic
            sensitivity_mat = [mu 0.0;0.0 -mu]
        elseif method == :explicit_staggered_grid_new
            sensitivity_mat = [mu 0.0;0.0 -mu]
        else
            sensitivity_mat = ForwardDiff.jacobian(end_flow_func, [rho_from, rho_to])
        end
        # @show end_flows, sensitivity_mat
        _assemble_pipe_solve_results!(ts, fr_node, to_node, end_flows, sensitivity_mat, r, J)
    end

    return
end


function solve_pipe_state!(ts::TransientSimulator, method::Symbol, pipe_id::Int64,rho_from::T,rho_to::T,inertial_flag = zero(T))::Vector{T} where {T<:Real}


    if method == :implicit_parabolic
        end_flows = _solve_pipe_state_parabolic!(ts, pipe_id, rho_from, rho_to, inertial_flag)

    elseif method == :explicit_hyperbolic
        end_flows = _solve_pipe_state_maccormack!(ts, pipe_id, rho_from, rho_to, inertial_flag)
    elseif method == :explicit_staggered_grid_new
        end_flows =  _solve_pipe_state_staggered_grid!(ts,pipe_id,rho_from,rho_to,inertial_flag)
    end

    return end_flows
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
        # println(iter, ":", residual)
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


function form_nodal_pressure_vector(ts::TransientSimulator)::Vector{Float64}
    key_array = sort(collect(keys(ref(ts, :node))))
    pressure_vector = Vector{Float64}(undef, length(key_array))
    for i in eachindex(key_array)
        pressure_vector[i] = ref(ts, :node, key_array[i])["pressure"]
    end
    return pressure_vector
end

function solve_newton_basic!(
    x::Vector{T},
    residual_fun!::Function,
    Jacobian_fun!::Function;
    tol::Real = 1e-8,
    max_iter::Int = 20) where {T<:Real}
    n = length(x)
    residual = zeros(T, n)
    J = (T == Float64) ? spzeros(T, n, n) : zeros(T, n, n)

    for iter = 1:max_iter
        fill!(residual, zero(T))
        residual_fun!(residual, x)
        res_norm = maximum(abs, residual)

        if res_norm <= tol
            return x, true, iter, res_norm
        end

        if J isa SparseMatrixCSC
            fill!(J.nzval, zero(T))
        else
            fill!(J, zero(T))
        end
        Jacobian_fun!(J, x)
        delta_x = J \ (-residual)
        x .+= delta_x

        step_norm = maximum(abs, delta_x)
        if step_norm <= tol
            return x, true, iter, res_norm
        end
    end

    fill!(residual, zero(T))
    residual_fun!(residual, x)
    return x, false, max_iter, maximum(abs, residual)
end

function check_limits(ts::TransientSimulator, x_node::Vector{Float64})

    pressure_min = params(ts, :minimum_pressure_limit) / nominal_values(ts, :pressure)
    rho_min = (pressure_min > 0) ? get_density(ts, pressure_min) : 0
    pressure_max = params(ts, :maximum_pressure_limit) / nominal_values(ts, :pressure)
    rho_max = get_density(ts, pressure_max)
    delta_x_node = zeros(length(x_node))
    violated_nodes = Vector{Int64}()
    for (index, x)  in enumerate(x_node)
        if x < rho_min
            push!(violated_nodes, index)
            delta_x_node[index] = rho_min - x
        elseif x > rho_max
            push!(violated_nodes, index)
            delta_x_node[index] = rho_max - x
        end
    end
    
    # check if any non-flow control compressor has ends i, j for which both i, j have delta_x_node nonzero.
    # if so, then no load adjustment possible because this the compressor equation cannot be satisfied unless it is the unlikely instance where these values are compatible with compressor equation.

    # if only one end of compressor, say i,  has an imposed density change, the compressor eqn now imposes a change at the other end  j too. this change  can be computed, but this means that at j, the new density after change must also be within limits. if not, then load adjustment is not possible. if yes, then load adjustment is possible with change at i and j. (this is step where this method could go wrong)

    # the systematic way to do this is to set violated nodes as new slack nodes  with limiting densities. solve problem agan. if no violations, compute these new "slack injections" at violated nodes
        
    #  check if residual of compressor eqns for x_node_new = x_node + delta_x_node is non-zero. is achievable by load adjustment at violated nodes. This requires evaluating J * delta_x_node where J is the Jacobian of the system w.r.t nodal densities, and checking if the required change in injection at violated nodes is compatible with compressor constraints. If not compatible, then throw error that load adjustment cannot fix density violation. If compatible, then update nodal injections at violated nodes accordingly to fix density violation. Note that this load adjustment step is not guaranteed to work and is a heuristic to try to fix density violations when they occur.
    # now evaluate J * delta_x_node

    if !isempty(violated_nodes)
        throw(DomainError("Density bound ($rho_min, $rho_max) violated at following node(s) $violated_nodes")
        )
    end
    return
end

# problem with load reduction is that it is not local and does not guarantee success.
# our aim is to reset nodal densities to the limit and check if this can be achieved by adjusting injections at violated nodes. how so ?
# suppose n nodes exceed bounds..WLOG lower bound. If you set them all to rho_min, then in the discrete system you solved, 
# you want equations at these nodes to be satisfied  with compensatory effect of change in nodal injection (load reduction).
# That is, you evaluate J delta_rho in general.  Note that almost all eqns are linear, so most rows of J sre const. For the linear eqns, J delta_rho represents the required change in injection for given changes in nodal density. JHowever, the problem here is that if some of these changes happen at nodes incident by a compressor, then you must have 
# p(rho_i + delta_rho_i) = alpha * p(rho_j + delta_rho_j). If you are given the increments at both ends, and  proposed density changes are compatible with compressor then load reduction is possible. else not. Alternatively, density changes at one end of the compressor imply a change at other end too. 