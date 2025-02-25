
using GasTranSim
using JSON

""" 
This function converts a steady state control to a ramp boundary condition 
to perform a transient simulation. The transient simulation will converge 
to the steady state given enough time. 
"""
function create_datafile_for_steady_state_time_marching(folder_in::AbstractString, folder_out::AbstractString; t_initial::Float64=0.0, t_ramp::Float64=21600.0, t_final::Float64=86400.0) 

    # folder_in = "/Users/shrirams/Documents/GasSteadySim.jl/examples/data/GasLib-40/"
    # folder_out = "./GasLib-40/"

    network_file = folder_in * "network_new.json"
    bc_steady_file = folder_in * "bc_new.json"
    params_steady_file = folder_in * "params_new.json"
    
    t_initial = 0.0
    t_ramp = 3600.0 * 6.0
    t_final = 3600 * 24.0 
    p_init = 1e5

    ic_filename = folder_out * "ic_ramp.json"
    bc_filename =  folder_out * "bc_ramp.json"
    params_filename = folder_out * "params_ramp.json"

    bc_ramp = Dict{String, Any}(
        "boundary_pslack" => Dict(), 
        "boundary_nonslack_flow" => Dict(), 
        "boundary_compressor" => Dict()
    )

    bc_steady_data = parse_json(bc_steady_file)

    for (_, p_val) in get(bc_steady_data, "boundary_pslack", [])
         p_init = p_val
        break
    end

    for (key, p_slack) in get(bc_steady_data, "boundary_pslack", [])
        bc_ramp["boundary_pslack"][key] = Dict{String, Any}()
        bc_ramp["boundary_pslack"][key]["time"] = [t_initial, t_ramp, t_final] 
        bc_ramp["boundary_pslack"][key]["value"] = [p_init, p_slack, p_slack]
    end

    for (key, flow_val) in get(bc_steady_data, "boundary_nonslack_flow", [])
        bc_ramp["boundary_nonslack_flow"][key] = Dict{String, Any}()
        bc_ramp["boundary_nonslack_flow"][key]["time"] = [t_initial, t_ramp, t_final] 
        bc_ramp["boundary_nonslack_flow"][key]["value"] = [0, flow_val, flow_val]
    end

    for (key, compressor_dict) in get(bc_steady_data, "boundary_compressor", [])
        bc_ramp["boundary_compressor"][key] = Dict{String, Any}()
        bc_ramp["boundary_compressor"][key]["time"] = [t_initial, t_ramp, t_final]
        ctrl = compressor_dict["control_type"]
        val = compressor_dict["value"]
        bc_ramp["boundary_compressor"][key]["control_type"] = [ctrl, ctrl, ctrl]
        bc_ramp["boundary_compressor"][key]["value"] = [1.0, val, val]
    end


    open(bc_filename, "w") do f 
        JSON.print(f, bc_ramp, 2)
    end

    network_data  = parse_json(network_file)

    ic_data  = Dict{String, Any}(
        "nodal_pressure" => Dict(), 
        "pipe_flow" => Dict(), 
        "compressor_flow" => Dict())


    for (node_id, _) in network_data["nodes"]
        ic_data["nodal_pressure"][node_id] = p_init
    end

    for (pipe_id, _) in network_data["pipes"]
        ic_data["pipe_flow"][pipe_id] = 0
    end

    for (comp_id, _) in network_data["compressors"]
        ic_data["compressor_flow"][comp_id] = 0
    end

    open(ic_filename, "w") do f
        JSON.print(f, ic_data, 2)
    end

    params_steady_data = parse_json(params_steady_file)

    params_dict = Dict{String, Any}("simulation_params" => Dict())

    params_dict["simulation_params"] = params_steady_data["params"]
    params_dict["simulation_params"]["Initial time"] = t_initial
    params_dict["simulation_params"]["Final time"] = t_final
    params_dict["simulation_params"]["Discretization time step"] = (t_final - t_initial)/1000
    params_dict["simulation_params"]["Output dx"] =  1000
    params_dict["simulation_params"]["Courant number (must be between 0 and 1, recommended value is 0.9"] = 0.5
    params_dict["simulation_params"]["Save final state"]=  1
    params_dict["simulation_params"]["Output dt"] =  3600

    open(params_filename, "w") do f
        JSON.print(f, params_dict, 2)
    end

    
end 

""" 
This function reads  a bc.json data file for a transient simulation and creates a steady-state boundary condition file based on the value at the first time point. 
"""
function create_datafile_for_steady_state_simulation_from_transient_datafile(folder::AbstractString)

    bc_steady = Dict{String, Any}(
    "boundary_pslack" => Dict(), 
    "boundary_nonslack_flow" => Dict(), 
    "boundary_compressor" => Dict())

    bc_file =  folder * "bc.json"
    bc_data = parse_json(bc_file)

    for (key, slack_dict) in bc_data["boundary_pslack"]
        bc_steady["boundary_pslack"][key] = slack_dict["value"][1]
    end
    for (key, flow_dict) in bc_data["boundary_nonslack_flow"]
        bc_steady["boundary_nonslack_flow"][key] = flow_dict["value"][1]
    end

    for (key, compressor_dict) in bc_data["boundary_compressor"]
        bc_steady["boundary_compressor"][key] = Dict{String, Any}()
        bc_steady["boundary_compressor"][key]["control_type"] = compressor_dict["control_type"][1]
        bc_steady["boundary_compressor"][key]["value"] = compressor_dict["value"][1]
    end


    bc_steady_file = folder * "bc_steady.json"
    open(bc_steady_file, "w") do f 
        JSON.print(f, bc_steady, 2)
    end

    params_file =  folder * "params.json"
    params_data = parse_json(params_file)

    params_dict = Dict{String, Any}("simulation_params" => Dict())
    params_dict["simulation_params"] = Dict("Specific heat capacity ratio"=> 1.4,
    "Temperature (K):"=> 288.70599999999996,
    "Gas specific gravity (G):"=> 0.6, 
    "units (SI = 0, standard = 1)" => 0)
    if haskey(params_data["simulation_params"], "nominal_pressure")
        params_dict["simulation_params"]["nominal_pressure"] = params_data["simulation_params"]["nominal_pressure"]
    end
    if haskey(params_data["simulation_params"], "nominal_velocity")
        params_dict["simulation_params"]["nominal_velocity"] = params_data["simulation_params"]["nominal_velocity"]
    end
    if haskey(params_data["simulation_params"], "nominal_length")
        params_dict["simulation_params"]["nominal_length"] = params_data["simulation_params"]["nominal_length"]
    end
    if haskey(params_data["simulation_params"], "nominal_density")
        params_dict["simulation_params"]["nominal_density"] = params_data["simulation_params"]["nominal_density"]
    end
    params_steady_file = folder * "params_steady.json"
    open(params_steady_file, "w") do f 
        JSON.print(f, params_dict, 2)
    end


    return
end