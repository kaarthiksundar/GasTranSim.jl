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
    bc = build_bc(data)
    pu_pressure_to_pu_density, pu_density_to_pu_pressure = get_eos(nominal_values, params, eos)

    ts = TransientSimulator(data,
        ref,
        initialize_solution(data, params),
        nominal_values,
        params,
        bc,
        pu_pressure_to_pu_density,
        pu_density_to_pu_pressure
    )

    add_pipe_grid_to_ref!(ts)
    add_node_level_flag!(ts)
    initialize_nodal_state!(ts)
    initialize_pipe_state!(ts)
    return ts
end


function add_pipe_grid_to_ref!(ts::TransientSimulator)
    for (key, pipe) in ref(ts, :pipe)

        # CFL condition a*dt/dx <= 1 => dx >= a*dt
        num_segments = pipe["length"] / (1.0 * params(ts, :dt)) 
        
        if num_segments < 1
            @error "Time step too large for CFL condition. Spatial discretization failed"
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
        ref(ts, :node, node_id)["is_level_2"] = false
        return
    end
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
    end
    return
end

function initialize_pipe_state!(ts::TransientSimulator)
    for (key, pipe) in ref(ts, :pipe)
        fill!(pipe["mass_flux_profile"], pipe["initial_mass_flux"])
        fr_node = pipe["fr_node"]
        to_node = pipe["to_node"]
        initial_fr_pressure = ref(ts, :node, fr_node, "pressure")
        initial_to_pressure = ref(ts, :node, to_node, "pressure")
        density_at_first_sq = get_density(ts, initial_fr_pressure) ^ 2
        density_at_last_sq = get_density(ts, initial_to_pressure) ^ 2
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
