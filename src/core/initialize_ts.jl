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
    method::Symbol = :explicit_staggered_grid,
    # method::Symbol = :implicit_parabolic,
)::TransientSimulator
    validate_method_contract!(method)
    params, nominal_values = process_data!(data)
    params[:method] = method
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

    initialize_pipe_grid!(ts, method)
    add_eqn_number_for_nodes!(ts)
    initialize_nodal_state!(ts)
    initialize_pipe_state!(ts, method)
    initialize_compressor_state!(ts)
    test_ic(ts)
    return ts
end
function initialize_nodal_state!(ts::TransientSimulator)
    for (key, _) in ref(ts, :node)
        pressure = initial_nodal_pressure(ts, key)
        ref(ts, :node, key)["pressure"] = pressure
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
    solve_compressor_flows!(ts, lin_system)
    return
end

function test_ic(ts::TransientSimulator)
    err_node = 0.0
    # error in nodal flow balance conditions
    for (node_id, _) in ref(ts, :node)
        if ref(ts, :node, node_id)["is_slack"] == 1
            continue
        end
        _, withdrawal = control(ts, :node, node_id, 0)

        flow_p = add_up_pipe_flows_at_node(node_id, ts)
        flow_c = add_up_compressor_flows_at_node(node_id, ts)
        term = flow_p + flow_c - withdrawal # withdrawal positive in input data already
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
    skip_nodes::Set{Int64})::Vector{Int64}
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
    slack_nodes::Vector{Int64})::Vector{Int64}
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


function form_matrix_for_compressor_flow_solve(ts::TransientSimulator)::Union{Tuple{SparseMatrixCSC{Float64,Int64},Vector{Int64},Vector{Int64}},Nothing}
    @info "Assuming that compressors alone do not form cycles. Else compressor flows cannot be determined uniquely"
    compressors = get(ref(ts), :compressor, Dict())
    isempty(compressors) && return nothing

    compressor_ids = sort(collect(keys(compressors)))
    num_compressors = length(compressor_ids)

    # Each compressor needs one unique non-slack endpoint equation; build candidates first.
    candidate_nodes = Dict{Int64,Vector{Int64}}()
    for c_id in compressor_ids
        fr_node = ref(ts, :compressor, c_id, "fr_node")
        to_node = ref(ts, :compressor, c_id, "to_node")

        cands = Int64[]
        (ref(ts, :node, fr_node)["is_slack"] == 0) && push!(cands, fr_node)
        if (ref(ts, :node, to_node)["is_slack"] == 0) && (to_node != fr_node)
            push!(cands, to_node)
        end

        isempty(cands) &&
            throw(NetworkException("Cannot solve for flow in Compressor $c_id: both endpoints are slack"))

        candidate_nodes[c_id] = cands
    end

    # Order-independent assignment via backtracking on compressors with fewer choices first.
    assignment_order = sort(compressor_ids, by = c_id -> length(candidate_nodes[c_id]))
    assigned_node = Dict{Int64,Int64}()
    used_nodes = Set{Int64}()

    function _assign_rhs_node!(k::Int64)::Bool
        k > length(assignment_order) && return true
        c_id = assignment_order[k]

        for node_id in candidate_nodes[c_id]
            node_id in used_nodes && continue
            assigned_node[c_id] = node_id
            push!(used_nodes, node_id)

            _assign_rhs_node!(k + 1) && return true

            delete!(assigned_node, c_id)
            delete!(used_nodes, node_id)
        end

        return false
    end

    _assign_rhs_node!(1) ||
        throw(NetworkException("Cannot solve compressor flows: no order-independent unique non-slack node assignment exists"))

    rhs_index = [assigned_node[c_id] for c_id in compressor_ids]

    A = spzeros(num_compressors, num_compressors)
    comp_col = Dict{Int64,Int64}(c_id => j for (j, c_id) in enumerate(compressor_ids))

    for i = 1:length(rhs_index)
        node_id = rhs_index[i]

        for c_id in ref(ts, :incoming_compressors, node_id)
            haskey(comp_col, c_id) || continue
            A[i, comp_col[c_id]] = 1.0
        end
        for c_id in ref(ts, :outgoing_compressors, node_id)
            haskey(comp_col, c_id) || continue
            A[i, comp_col[c_id]] = -1.0
        end
    end

    return A, rhs_index, compressor_ids
end

function solve_compressor_flows!(ts::TransientSimulator, lin_system::Union{Tuple,Nothing})
    if isnothing(lin_system)
        return
    else
        A, rhs_index, compressor_ids = lin_system
    end
    b = zeros(length(rhs_index))

    for i = 1:length(rhs_index)
        node_id = rhs_index[i]
        ctrl_type, ctrl_val = control(ts, :node, node_id, ref(ts, :current_time))
        flow = add_up_pipe_flows_at_node(node_id, ts) 
        b[i] = -flow + ctrl_val # withdrawal positive in input data already
    end
    x = A \ b

    for i = 1:length(compressor_ids)
        c_id = compressor_ids[i]
        ref(ts, :compressor, c_id)["flow"] = x[i]
    end

    return
end

function calculate_slack_injections!(ts::TransientSimulator)
    for (node_id, node) in ref(ts, :node)
        if node["is_slack"] == 1
            
            flow_p = add_up_pipe_flows_at_node(node_id, ts)
            flow_c = add_up_compressor_flows_at_node(node_id, ts)
            ref(ts, :node, node_id)["injection"] = -(flow_p + flow_c)
        end
    end
    return
end


function add_up_pipe_flows_at_node(node_id::Int64,
    ts::TransientSimulator)::Real
    out_p = ref(ts, :outgoing_pipes, node_id)
    in_p = ref(ts, :incoming_pipes, node_id)
    flow = 0.0 
    pipe = ref(ts, :pipe)

    for p in out_p
        flow -= pipe[p]["fr_mass_flux"] * pipe[p]["area"]
    end
    for p in in_p
        flow += pipe[p]["to_mass_flux"] * pipe[p]["area"]
    end
    return flow
end

function add_up_compressor_flows_at_node(node_id::Int64,
    ts::TransientSimulator)::Real
    out_c = ref(ts, :outgoing_compressors, node_id)
    in_c = ref(ts, :incoming_compressors, node_id)
    flow = 0.0 
    t = ref(ts, :current_time)
    for ci in in_c
        ctr, cmpr_val = control(ts, :compressor, ci, t)
        if ctr == flow_control
            flow += cmpr_val # inflow positive
        end
    end
    for co in out_c
        ctr, cmpr_val = control(ts, :compressor, co, t)
        if ctr == flow_control
            flow += (-1.0 * cmpr_val) # outflow negative
        end
    end
    return flow
end
