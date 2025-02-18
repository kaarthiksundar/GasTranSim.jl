function get_data(data_folder::AbstractString; 
    case_name::AbstractString="", 
    case_types::Vector{Symbol}=Symbol[])::Dict{String,Any}
    return parse_data(data_folder; case_name=case_name, case_types=case_types)
end

function initialize_simulator(data_folder::AbstractString;
    case_name::AbstractString="", 
    case_types::Vector{Symbol}=Symbol[],
    kwargs...)::TransientSimulator
    data = parse_data(data_folder; case_name=case_name, case_types=case_types)
    return initialize_simulator(data; kwargs...)
end

function initialize_simulator(data::Dict{String,Any}; eos::Symbol=:ideal)::TransientSimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    ref = build_ref(data, ref_extensions= [
        add_pipe_info_at_nodes!,
        add_compressor_info_at_nodes!,
        add_current_time_to_ref!]
    )
    ref[:current_time] = params[:t_0]

    ts = TransientSimulator(data,
        ref,
        initialize_solution(data),
        nominal_values,
        params,
        build_ic(data),
        build_bc(data),
        get_eos(eos)...
    )

    add_pipe_grid_to_ref!(ts)
    add_node_level_flag!(ts)
    add_node_ordering_for_compressor_flow_computation!(ts)
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
        num_segments = c_inv * (pipe["length"] * params(ts, :courant_number)) / params(ts, :dt) 
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

function _evaluate_level_of_node!(ts::TransientSimulator, node_id::Int64)
    if !haskey(ref(ts), :compressor)
        ref(ts, :node, node_id)["level"] = 0 
        return
    end

    # forbid topology if compressor delivers to slack node
    if ref(ts, :node, node_id)["is_slack"] == true && length(ref(ts, :incoming_compressors, node_id)) > 0
        throw(NetworkException("Compressor delivering to slack node"))
    end

    if length(ref(ts, :incoming_compressors, node_id)) > 1
        num_slack_fr_nodes = 0
        for ci in ref(ts, :incoming_compressors, node_id)
            node_across_ci = ref(ts, :compressor, ci, "fr_node")
            if ref(ts, :node, node_across_ci)["is_slack"] == true
                num_slack_fr_nodes += 1
            end
        end
        if num_slack_fr_nodes > 1
            throw(NetworkException("Two slack nodes connected by two compressors in series"))
        end
    end

    if length(ref(ts, :incoming_compressors, node_id)) + length(ref(ts, :outgoing_compressors, node_id)) == 0
        ref(ts, :node, node_id)["level"] = 0 
        return
    end

    if length(ref(ts, :incoming_compressors, node_id)) + length(ref(ts, :outgoing_compressors, node_id)) == 1
        if length(ref(ts, :incoming_compressors, node_id)) == 1
            ref(ts, :node, node_id)["level"] = -1
        end
        if length(ref(ts, :outgoing_compressors, node_id)) == 1
            ref(ts, :node, node_id)["level"] =  1
        end
        return
    end

    for ci in ref(ts, :incoming_compressors, node_id)
        node_across_ci = ref(ts, :compressor, ci, "fr_node")
        if length(ref(ts, :incoming_compressors, node_across_ci)) + length(ref(ts, :outgoing_compressors, node_across_ci))  >  1
            throw(NetworkException("3 compressors are in series"))
        end
    end

    for ci in ref(ts, :outgoing_compressors, node_id)
        node_across_ci = ref(ts, :compressor, ci, "to_node")
        if length(ref(ts, :incoming_compressors, node_across_ci)) + length(ref(ts, :outgoing_compressors, node_across_ci))  > 1
            throw(NetworkException("3 compressors are in series"))
        end
    end

    ref(ts, :node, node_id)["level"] = 2
    return
end

function add_node_level_flag!(ts::TransientSimulator)
    for (node_id, _) in ref(ts, :node)
        _evaluate_level_of_node!(ts, node_id)
    end
    return
end

function add_node_ordering_for_compressor_flow_computation!(ts::TransientSimulator)
    compressors = get(ref(ts), :compressor, Dict())
    (isempty(compressors)) && (return)
    compressor_ids = keys(compressors) |> collect
    num_compressors = length(compressors)
    compressors_accounted_for = Vector{Int64}()
    compressors_remaining = Vector{Int64}()
    used_node_ids = Vector{Int64}()

    for c_id in compressor_ids
        fr_id = ref(ts, :compressor, c_id)["fr_node"]
        to_id = ref(ts, :compressor, c_id)["to_node"]
        if abs(ref(ts, :node, fr_id)["level"]) == 1 && ref(ts, :node, fr_id)["is_slack"] == false
            push!(compressors_accounted_for, c_id)
            push!(used_node_ids, fr_id)
            continue
        elseif abs(ref(ts, :node, to_id)["level"]) == 1 # to_id cannot be slack
            push!(compressors_accounted_for, c_id)
            push!(used_node_ids, to_id)
            continue
        else
            push!(compressors_remaining, c_id)
        end
    end

    if length(compressors_accounted_for) < num_compressors
        @debug "The flows in some compressors $compressors_remaining cannot be calculated without knowing  other compressor flows"
    end 

    ref(ts)[:node_ordering_for_compressor_flow_calculation] = used_node_ids 
    ref(ts)[:compressor_ids_second_round_calculation] = compressors_remaining
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
                    density_at_last_sq * (i - 1) * dL +
                    density_at_first_sq * (n - i) * dL
                ) for i in 1:n
            ]
            pipe["fr_minus_mass_flux"] = initial_mass_flux # (dx/2)
            pipe["to_minus_mass_flux"] = initial_mass_flux # L-(dx/2)
            pipe["fr_mass_flux"] = initial_mass_flux # -(dx/2)
            pipe["to_mass_flux"] = initial_mass_flux # L+(dx/2)
        else 
            flow_spl = initial_pipe_mass_flow(ts, key)
            pressure_spl = initial_pipe_pressure(ts, key)
            x_rho = LinRange(0, L, n)
            x_mid = x_rho[1:n-1] .+ dx/2.0
            get_coeffs(flow_spl)[1]
            pipe["mass_flux_profile"] = [
                get_coeffs(flow_spl)[1], 
                [flow_spl(x) for x in x_mid]..., 
                get_coeffs(flow_spl)[end]
            ] ./ area
            pipe["fr_minus_mass_flux"] = pipe["mass_flux_profile"][2]
            pipe["to_minus_mass_flux"] = pipe["mass_flux_profile"][end-1]
            pipe["fr_mass_flux"] = pipe["mass_flux_profile"][1]
            pipe["to_mass_flux"] = pipe["mass_flux_profile"][end]
            pipe["density_profile"][1:n] = [
                get_density(ts, pressure_spl(x)) for x in x_rho
            ]
        end 
    end
    return
end

function _compute_compressor_flows!(ts::TransientSimulator)

    for node_id in get(ref(ts), :node_ordering_for_compressor_flow_calculation, [])
        t = ref(ts, :current_time)
        ctrl_type, withdrawal = control(ts, :node, node_id, t)
        @assert ctrl_type == flow_control
        _, net_injection = _assemble_pipe_contributions_to_node(node_id, withdrawal, 1.0, ts)

        if ref(ts, :node, node_id)["level"] == 1
            # find the outgoing compressor
            out_c = ref(ts, :outgoing_compressors)[node_id]
            ref(ts, :compressor, out_c[1])["flow"] = net_injection
        end
        if ref(ts, :node, node_id)["level"] == -1
            # find the incoming compressor
            in_c = ref(ts, :incoming_compressors)[node_id]
            ref(ts, :compressor, in_c[1])["flow"] = -net_injection
        end
    end

    # Network topology conditions imply following calculation cannot fail
    # Flows of other compressors at node must be known
    # Compressor has unknown flow must mean one end is slack and other end is level 2
    # Both ends cannot be level 2 or slack, if one end was level 1 non-slack, would already be finished
    for c_id in get(ref(ts), :compressor_ids_second_round_calculation, [])
        to_id = ref(ts, :compressor, c_id)["to_node"]
        # to_id cannot be slack, so focus on that
        t = ref(ts, :current_time)
        ctrl_type, withdrawal = control(ts, :node, to_id, t)
        @assert ctrl_type == flow_control
        _, net_injection = _assemble_pipe_contributions_to_node(to_id, withdrawal, 1.0, ts)
        in_c = ref(ts, :incoming_compressors)[to_id]
        out_c = ref(ts, :outgoing_compressors)[to_id]
        for id in  in_c 
            if id == c_id # c_id is incoming
                continue
            end
            net_injection += ref(ts, :compressor, id)["flow"]
        end
        for id in  out_c
            net_injection -= ref(ts, :compressor, id)["flow"]
        end
        ref(ts, :compressor, c_id)["flow"] = -net_injection
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
    _compute_compressor_flows!(ts)
    return
end

function test_ic(ts::TransientSimulator)
    err_node = 0.0
    for (node_id, junction) in ref(ts, :node)
        if ref(ts, :node, node_id)["is_slack"] == 1
            continue
        end
        ctrl_type, term = control(ts, :node, node_id, 0)

        out_p = ref(ts, :outgoing_pipes, node_id)
        in_p = ref(ts, :incoming_pipes, node_id)
        for i in out_p
            term += ref(ts, :pipe, i, "fr_mass_flux") * ref(ts, :pipe, i, "area") #withdrawal positive
        end

        for i in in_p
            term -= ref(ts, :pipe, i, "fr_mass_flux") * ref(ts, :pipe, i, "area") #withdrawal positive
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
    for (key, compressor) in ref(ts, :compressor)
        ctrl_type, alpha = control(ts, :compressor, key, 10)
        if  ctrl_type == 0
            fr_pr = initial_nodal_pressure(ts, ref(ts, :compressor, key, "fr_node"))
            to_pr = initial_nodal_pressure(ts, ref(ts, :compressor, key, "to_node"))
            term = abs(alpha * fr_pr - to_pr)
            err_c = max(err_c, term)
        end
    end


    err_p = 0.0
    for (key, pipe) in ref(ts, :pipe)
        fr_pr = initial_nodal_pressure(ts, ref(ts, :pipe, key, "fr_node"))
        to_pr = initial_nodal_pressure(ts, ref(ts, :pipe, key, "to_node"))
        pr_mean = (fr_pr + to_pr) / 2
        rho_mean = get_density(ts, pr_mean)
        term = (fr_pr  - to_pr) * rho_mean - pipe["length"] *  pipe["friction_factor"] / (2 * pipe["diameter"] ) * pipe["fr_mass_flux"] * abs(pipe["fr_mass_flux"])
        err_p = max(err_p, abs(term))
    end
    err = max(err_node, err_c, err_p)
    @info "Error in initial condition is $err"
    
    return
end