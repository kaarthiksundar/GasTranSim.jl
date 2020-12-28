function initialize_simulator(file::AbstractString)::TransientSimulator
    data = parse_json(file)
    return initialize_simulator(data)
end 

function initialize_simulator(data::Dict{String,Any})::TransientSimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    ref = build_ref(data, ref_extensions= [add_pipe_info_at_nodes!, add_compressor_info_at_nodes!,
        add_current_time_to_ref!]) 
    ref[:current_time] = params[:t_0]
    bc = build_bc(data)
    pu_pressure_to_pu_density = x -> x 
    pu_density_to_pu_pressure = x -> x 
    pu_dp_drho = x-> x^0

    return TransientSimulator(data, 
        ref, 
        Dict{String,Any}(), 
        nominal_values, 
        params, 
        bc, 
        pu_pressure_to_pu_density, 
        pu_density_to_pu_pressure,
        pu_dp_drho
        ) 
end 



function add_grid_to_ref!(ts::TransientSimulator)
    for (key, pipe) in ts.ref[:pipe]
        n = ceil(Int64, ts.ref[:pipe][key]["length"] / (1 * ts.params[:dt]) ) 
        ts.ref[:pipe][key]["num_discretization_points"] = n
        ts.ref[:pipe][key]["dx"] = ts.ref[:pipe][key]["length"] / (n-1)
        ts.ref[:pipe][key]["density_profile"] = zeros(Float64, n)
        ts.ref[:pipe][key]["mass_flux_profile"] = zeros(Float64, n+1)
    end
    return
end

function initialize_vertex!(ts::TransientSimulator)

    for (key, junction) in ts.ref[:node]
        ts.ref[:node][key]["pressure"] = ts.ref[:node][key]["initial_pressure"]
    end
    return
end

function initialize_pipes!(ts::TransientSimulator)
    for (key, pipe) in ts.ref[:pipe]
        fill!(pipe["mass_flux_profile"], pipe["initial_mass_flux"])
        r1_sq = get_density(ts, pipe["initial_fr_pressure"]) ^ 2 
        rn_sq = get_density(ts, pipe["initial_to_pressure"]) ^ 2
        n = pipe["num_discretization_points"]
        dx = pipe["dx"]
        L = pipe["length"]
        for i = 1: n
            pipe["density_profile"][i] = sqrt( rn_sq * (i-1)* dx / L  + r1_sq *  (n - i) * dx / L )  
        end
        pipe["fr_minus_mass_flux"] = pipe["initial_mass_flux"]
        pipe["to_minus_mass_flux"] = pipe["initial_mass_flux"]
        pipe["fr_mass_flux"] = pipe["initial_mass_flux"]
        pipe["to_mass_flux"] = pipe["initial_mass_flux"]
    end
    return
end


# 3 conditions that network must satisfy
# 1. if we visit each node of level 0/1 (assuming a star configuration) and update immediate neighbours, 
 # all vertices should be eventually covered in finite steps. 
# 2. Node A cannot be delivery end of two  compressors with delivery pressure controls
# 3. node A cannot be both  slack node and delivery end of  compressor
