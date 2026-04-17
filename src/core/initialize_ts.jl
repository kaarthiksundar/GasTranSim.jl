function get_data(
    data_folder::AbstractString;
    case_name::AbstractString = "",
    case_types::Vector{Symbol} = Symbol[],
)::Dict{String,Any}
    return parse_data(data_folder; case_name = case_name, case_types = case_types)
end

function initialize_simulator(
    data_folder::AbstractString;
    case_name::AbstractString = "",
    case_types::Vector{Symbol} = Symbol[],
    kwargs...,
)::TransientSimulator
    data = parse_data(data_folder; case_name = case_name, case_types = case_types)
    return initialize_simulator(data; kwargs...)
end

function initialize_simulator(
    data::Dict{String,Any};
    eos::Symbol = :ideal,
)::TransientSimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    ref = build_ref(
        data,
        ref_extensions = [
            add_pipe_info_at_nodes!,
            add_compressor_info_at_nodes!,
            add_current_time_to_ref!,
        ],
    )
    ref[:current_time] = params[:t_0]

    ts = TransientSimulator(
        data,
        ref,
        initialize_solution(data),
        nominal_values,
        params,
        build_ic(data),
        build_bc(data),
        get_eos(eos)...,
    )

    add_pipe_grid_to_ref!(ts)
    add_eqn_number_for_nodes!(ts)
    initialize_nodal_state!(ts)
    initialize_pipe_state!(ts)
    initialize_compressor_state!(ts)
    test_ic(ts)
    return ts
end


function add_pipe_grid_to_ref!(ts::TransientSimulator)
    for (key, pipe) in ref(ts, :pipe)
        # CFL condition c*dt/dx <= 0.9 => dx >= c*dt/0.9
        # with nondim dt, dx, we have nondim_dt/ nondim_dx < = 0.9 * mach_no
        c_inv = nominal_values(ts, :mach_num)
        num_segments =
            c_inv * (pipe["length"] * params(ts, :courant_number)) / params(ts, :base_dt)
        n = 1
        if num_segments < 1
            throw(CFLException(string(key)))
        else
            n = floor(Int64, num_segments) + 1
        end
        ref(ts, :pipe, key)["num_discretization_points"] = n
        ref(ts, :pipe, key)["dx"] = pipe["length"] / (n-1)
        ref(ts, :pipe, key)["density_profile"] = zeros(Float64, n)
        ref(ts, :pipe, key)["mass_flux_profile"] = zeros(Float64, n+1)
    end
    return
end

function initialize_nodal_state!(ts::TransientSimulator)
    for (key, _) in ref(ts, :node)
        pressure = initial_nodal_pressure(ts, key)
        ref(ts, :node, key)["pressure"] = pressure
    end
    return
end


function initialize_pipe_state!(ts::TransientSimulator)
    is_steady = false
    # we are assuming initial pipe pressures will not be provided only for steady initial conditions 
    # the case where initial pipe flow is unsteady, explicitly computing initial pipe pressures is very tedious and hence not done 
    # when initial pipe flow in unsteady, we assume initial pipe pressure are provided.
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
            density_at_first_sq = get_density(ts, initial_fr_pressure) ^ 2
            density_at_last_sq = get_density(ts, initial_to_pressure) ^ 2
            dL = dx / L
            pipe["density_profile"][1:n] = [
                sqrt(
                    density_at_last_sq * (i - 1) * dL + density_at_first_sq * (n - i) * dL,
                ) for i = 1:n
            ]
            pipe["fr_minus_mass_flux"] = initial_mass_flux # (dx/2)
            pipe["to_minus_mass_flux"] = initial_mass_flux # L-(dx/2)
            pipe["fr_mass_flux"] = initial_mass_flux # -(dx/2)
            pipe["to_mass_flux"] = initial_mass_flux # L+(dx/2)
        else
            flow_spl = initial_pipe_mass_flow(ts, key)
            pressure_spl = initial_pipe_pressure(ts, key)
            x_rho = LinRange(0, L, n)
            x_mid = x_rho[1:(n-1)] .+ dx/2.0
            get_coeffs(flow_spl)[1]
            pipe["mass_flux_profile"] =
                [
                    get_coeffs(flow_spl)[1],
                    [flow_spl(x) for x in x_mid]...,
                    get_coeffs(flow_spl)[end],
                ] ./ area
            pipe["fr_minus_mass_flux"] = pipe["mass_flux_profile"][2]
            pipe["to_minus_mass_flux"] = pipe["mass_flux_profile"][end-1]
            pipe["fr_mass_flux"] = pipe["mass_flux_profile"][1]
            pipe["to_mass_flux"] = pipe["mass_flux_profile"][end]
            pipe["density_profile"][1:n] = [get_density(ts, pressure_spl(x)) for x in x_rho]
        end
    end
    return
end

function _solve_compressor_flows!(ts::TransientSimulator, lin_system::Union{Tuple,Nothing})
    if isnothing(lin_system)
        return
    else
        A, rhs_index = lin_system
    end
    b = zeros(length(rhs_index))

    for i = 1:length(rhs_index)
        node_id = rhs_index[i]
        ctrl_type, ctrl_val = control(ts, :node, node_id, ref(ts, :current_time))
        _ , bi = _assemble_pipe_contributions_to_node_new(node_id, ctrl_val, ts)
        b[i] = -bi
    end
    x  = A \ b

    for c_id = 1:length(rhs_index)
        ref(ts, :compressor, c_id)["flow"] = x[c_id]
    end
    
    return
end


function initialize_compressor_state!(ts::TransientSimulator)
    isempty(get(ref(ts), :compressor, [])) && return
    if has_compressor_initial_flows(ts)
        for (key, _) in ref(ts, :compressor)
            flow = initial_compressor_flow(ts, key)
            ref(ts, :compressor, key)["flow"] = flow
        end
        return
    end
    @info "compressor does not have initial flows, computing them"
    lin_system = form_matrix_for_compressor_flow_solve(ts)
    _solve_compressor_flows!(ts, lin_system)
    return
end

function test_ic(ts::TransientSimulator)
    err_node = 0.0
    # error in nodal flow balance conditions
    for (node_id, _) in ref(ts, :node)
        if ref(ts, :node, node_id)["is_slack"] == 1
            continue
        end
        _, term = control(ts, :node, node_id, 0)

        out_p = ref(ts, :outgoing_pipes, node_id)
        in_p = ref(ts, :incoming_pipes, node_id)
        for i in out_p
            term += ref(ts, :pipe, i, "fr_mass_flux") * ref(ts, :pipe, i, "area") # withdrawal positive
        end

        for i in in_p
            term -= ref(ts, :pipe, i, "fr_mass_flux") * ref(ts, :pipe, i, "area") # withdrawal positive
        end

        out_c = ref(ts, :outgoing_compressors, node_id)
        in_c = ref(ts, :incoming_compressors, node_id)

        for i in out_c
            term += ref(ts, :compressor, i)["flow"] #withdrawal positive
        end

        for i in in_c
            term -= ref(ts, :compressor, i)["flow"]
        end

        err_node = max(err_node, abs(term))
    end

    err_c = 0.0
    # error in compressor ratio equation (only if ctrl_type == 0 or 1)
    for (key, _) in get(ref(ts), :compressor, [])
        ctrl_type, val = control(ts, :compressor, key, 0)
        term = 0.0
        if ctrl_type == 0
            fr_pr = initial_nodal_pressure(ts, ref(ts, :compressor, key, "fr_node"))
            to_pr = initial_nodal_pressure(ts, ref(ts, :compressor, key, "to_node"))
            term = abs(val * fr_pr - to_pr)
        elseif ctrl_type == 1
            to_pr = initial_nodal_pressure(ts, ref(ts, :compressor, key, "to_node"))
            term = abs(val - to_pr)
        end
        err_c = max(err_c, term)
    end

    # pipe equation error 
    err_p = 0.0
    for (key, pipe) in ref(ts, :pipe)
        fr_pr = initial_nodal_pressure(ts, ref(ts, :pipe, key, "fr_node"))
        to_pr = initial_nodal_pressure(ts, ref(ts, :pipe, key, "to_node"))
        pr_mean = (fr_pr + to_pr) / 2
        rho_mean = get_density(ts, pr_mean)
        term =
            (fr_pr - to_pr) * rho_mean -
            pipe["length"] * pipe["friction_factor"] / (2 * pipe["diameter"]) *
            pipe["fr_mass_flux"] *
            abs(pipe["fr_mass_flux"])
        err_p = max(err_p, abs(term))
    end
    err = max(err_node, err_c, err_p)
    @debug "If the initial condition was steady, then the error in initial condition is $err"

    return
end



function _collect_compressor_component_nodes(
    ts::TransientSimulator,
    start_node::Int64,
    skip_nodes::Set{Int64},
)::Vector{Int64}
    component = Int64[]
    queue = Int64[start_node]
    seen = Set{Int64}([start_node])
    head = 1

    # Breadth-first traversal over nodes connected by non-flow-control compressors.
    while head <= length(queue)
        node_id = queue[head]
        head += 1
        push!(component, node_id)

        in_c = ref(ts, :incoming_compressors, node_id)
        out_c = ref(ts, :outgoing_compressors, node_id)

        for c in vcat(in_c, out_c)
            ctrl_type, _ = control(ts, :compressor, c, 0)
            (ctrl_type == flow_control) && continue

            fr = ref(ts, :compressor, c, "fr_node")
            to = ref(ts, :compressor, c, "to_node")
            for nbr in (fr, to)
                if (nbr in skip_nodes) || (nbr in seen)
                    continue
                end
                push!(seen, nbr)
                push!(queue, nbr)
            end
        end
    end

    return component
end

function add_eqn_number_for_nodes!(ts::TransientSimulator)
    node_ids = sort(collect(keys(ref(ts, :node))))
    assignment_vec = falses(length(node_ids))
    compressors = get(ref(ts), :compressor, Dict())

    # Only non-flow-control compressors contribute compressor equations.
    num_non_flow_compressors = 0
    for (comp_id, _) in compressors
        ctrl_type, _ = control(ts, :compressor, comp_id, 0)
        (ctrl_type == flow_control) && continue
        num_non_flow_compressors += 1
    end

    eqn_init_num = num_non_flow_compressors
    ts.ref[:num_compressors] = eqn_init_num

    # Nodes reachable from slack nodes through non-flow-control compressors do not
    # get nodal equations in this reduced formulation.
    slack_node_list = [n for (n, nd) in ref(ts, :node) if nd["is_slack"] == 1]
    nodes_no_eq_list = collect_no_eqn_nodes_from_slacks(ts, slack_node_list)
    no_eq_set = Set(nodes_no_eq_list)

    for node_id in nodes_no_eq_list
        ref(ts, :node, node_id)["eqn_number"] = NaN
        assignment_vec[node_id] = true
    end

    # Every connected component through non-flow-control compressors shares one
    # nodal equation number.
    for node_id in node_ids
        assignment_vec[node_id] && continue
        component = _collect_compressor_component_nodes(ts, node_id, no_eq_set)
        eqn_init_num += 1
        for cid in component
            ref(ts, :node, cid)["eqn_number"] = eqn_init_num
            assignment_vec[cid] = true
        end
    end

    if any(assignment_vec .== false)
        throw(NetworkException("Not all nodes were assigned an equation number"))
    end
    return
end

"""
Starting from slack nodes, traverse through non-flow-control compressors and
collect reachable non-slack nodes into `no_eqn_list`.
"""
function collect_no_eqn_nodes_from_slacks(
    ts::TransientSimulator,
    slack_nodes::Vector{Int64},
)::Vector{Int64}
    explore_list = copy(slack_nodes)
    no_eqn_list = Int64[]

    seen = Set{Int64}(slack_nodes)
    head = 1

    while head <= length(explore_list)
        node_id = explore_list[head]
        head += 1

        in_c = ref(ts, :incoming_compressors, node_id)
        out_c = ref(ts, :outgoing_compressors, node_id)

        for c in vcat(in_c, out_c)
            # Ignore flow-control compressors when building no-equation regions.
            ctrl_type, _ = control(ts, :compressor, c, 0)
            (ctrl_type == flow_control) && continue

            fr = ref(ts, :compressor, c, "fr_node")
            to = ref(ts, :compressor, c, "to_node")

            for nbr in (fr, to)
                if !(nbr in seen)
                    push!(seen, nbr)
                    push!(explore_list, nbr)
                    push!(no_eqn_list, nbr)
                end
            end
        end
    end

    slack_set = Set(slack_nodes)
    return [n for n in no_eqn_list if !(n in slack_set)]
end


function form_matrix_for_compressor_flow_solve(ts::TransientSimulator)::Union{Tuple{SparseMatrixCSC{Float64,Int64}, Vector{Int64}}, Nothing}
    @info "Assuming that compressors alone do not form cycles. Else compressor flows cannot be determined uniquely"
    compressors = get(ref(ts), :compressor, Dict())
    (isempty(compressors)) && return nothing
    num_compressors = length(compressors)
    A = spzeros(num_compressors, num_compressors)
    rhs_index = Vector{Int64}()

    for i = 1 : num_compressors
        fr_node = ref(ts, :compressor, i, "fr_node")
        to_node = ref(ts, :compressor, i, "to_node")
        
        if ref(ts, :node, fr_node)["is_slack"] == 0 && !(fr_node in rhs_index)
            push!(rhs_index, fr_node)
        elseif ref(ts, :node, to_node)["is_slack"] == 0 && !(to_node in rhs_index)
            push!(rhs_index, to_node)
        else
            throw(NetworkException("Cannot solve for flow in Compressor $i "))
        end
    end

    for i = 1: length(rhs_index)
        for c_id in ref(ts, :incoming_compressors, rhs_index[i])
            A[i, c_id] = 1.0
        end
        for c_id in ref(ts, :outgoing_compressors, rhs_index[i])
            A[i, c_id] = -1.0
        end
    end

    return A, rhs_index
end

function calculate_slack_injections!(ts::TransientSimulator)
    for (node_id, node) in ref(ts, :node)
        if node["is_slack"] == 1
            
            _ , flow = _assemble_pipe_contributions_to_node_new(node_id, 0.0, ts)
    
            out_c = ref(ts, :outgoing_compressors, node_id)
            in_c = ref(ts, :incoming_compressors, node_id)

            for i in out_c
                flow -= ref(ts, :compressor, i)["flow"] #withdrawal positive
            end

            for i in in_c
                flow += ref(ts, :compressor, i)["flow"]
            end

            ref(ts, :node, node_id)["injection"] = -flow
        end

    end
    return
end