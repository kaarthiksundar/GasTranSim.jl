function initialize_simulator(file::AbstractString)::TransientSimulator
    data = parse_json(file)
    return initialize_simulator(data)
end

function initialize_simulator(data::Dict{String,Any})::TransientSimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    ref = build_ref(data, ref_extensions= [
        add_pipe_info_at_nodes!,
        add_compressor_info_at_nodes!,
        add_current_time_to_ref!]
    )
    ref[:current_time] = params[:t_0]
    bc = build_bc(data)
    pu_pressure_to_pu_density = x -> x
    pu_density_to_pu_pressure = x -> x
    pu_dpressure_to_pu_ddensity = x-> x^0

    ts = TransientSimulator(data,
        ref,
        initialize_solution(data, params),
        nominal_values,
        params,
        bc,
        pu_pressure_to_pu_density,
        pu_density_to_pu_pressure,
        pu_dpressure_to_pu_ddensity
    )

    add_pipe_grid_to_ref!(ts)
    add_node_level_flag!(ts)
    initialize_nodal_state!(ts)
    initialize_pipe_state!(ts)
    return ts
end


function add_pipe_grid_to_ref!(ts::TransientSimulator)
    for (key, pipe) in ref(ts, :pipe)
        #TODO: account for courant number and for non-ideal effects
        n = ceil(Int64, pipe["length"] / (1.0 * params(ts, :dt)) )
        ref(ts, :pipe, key)["num_discretization_points"] = n
        ref(ts, :pipe, key)["dx"] = pipe["length"] / (n-1)
        ref(ts, :pipe, key)["density_profile"] = zeros(Float64, n)
        ref(ts, :pipe, key)["mass_flux_profile"] = zeros(Float64, n+1)
    end
    return
end

function _evaluate_level_of_node!(ts::TransientSimulator, node_id::Int64)
    if length(ref(ts, :incoming_compressors, node_id)) + length(ref(ts, :outgoing_compressors, node_id)) == 0
        ref(ts, :node, node_id)["is_level_2"] = false
        return
    end
    for ci in ref(ts, :incoming_compressors, node_id)
        node_across_ci = ref(ts, :compressor, ci, "fr_node")
        if length(ref(ts, :incoming_compressors, node_across_ci)) + length(ref(ts, :outgoing_compressors, node_across_ci))  > 1
            ref(ts, :node, node_id)["is_level_2"] = true
            return
        end
    end
    for ci in ref(ts, :outgoing_compressors, node_id)
        node_across_ci = ref(ts, :compressor, ci, "to_node")
        if length(ref(ts, :incoming_compressors, node_across_ci)) + length(ref(ts, :outgoing_compressors, node_across_ci))  > 1
            ref(ts, :node, node_id)["is_level_2"] = true
            return
        end
    end
    ref(ts, :node, node_id)["is_level_2"] = false
    return
end

function add_node_level_flag!(ts::TransientSimulator)
    for (node_id, node) in ref(ts, :node)
        _evaluate_level_of_node!(ts, node_id)
    end
    return
end


function initialize_nodal_state!(ts::TransientSimulator)
    for (key, node) in ref(ts, :node)
        pressure = node["initial_pressure"]
        ref(ts, :node, key)["pressure"] = pressure
        ref(ts, :node, key)["density"] = get_density(ts, pressure)
    end
    return
end

function initialize_pipe_state!(ts::TransientSimulator)
    for (key, pipe) in ref(ts, :pipe)
        fill!(pipe["mass_flux_profile"], pipe["initial_mass_flux"])
        density_at_first_sq = get_density(ts, pipe["initial_fr_pressure"]) ^ 2
        density_at_last_sq = get_density(ts, pipe["initial_to_pressure"]) ^ 2
        n = pipe["num_discretization_points"]
        dx = pipe["dx"]
        L = pipe["length"]
        dL = dx / L
        pipe["density_profile"][1:n] = [
            sqrt(
                    density_at_last_sq * (i - 1) * dL +
                    density_at_first_sq * (n - i) * dL
            ) for i in 1:n
        ]
        pipe["fr_minus_mass_flux"] = pipe["initial_mass_flux"]
        pipe["to_minus_mass_flux"] = pipe["initial_mass_flux"]
        pipe["fr_mass_flux"] = pipe["initial_mass_flux"]
        pipe["to_mass_flux"] = pipe["initial_mass_flux"]
    end
    return
end



