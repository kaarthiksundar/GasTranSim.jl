
using GasTranSim
using JSON

""" 
This function converts a steady state control to a ramp boundary condition 
to perform a transient simulation. The transient simulation will converge 
to the steady state given enough time. 
"""
function convert_steady_solution_to_ramp_bc(folder::AbstractString) 
    network_file = folder * "network.json"
    bc_steady_file = folder * "bc_steady.json"
    params_file = folder * "params.json"
    
    t_initial = 0.0
    t_ramp = 3600.0 * 6.0
    t_final = 3600 * 24.0 
    p_init = 1e5

    ic_filename = folder * "ic_ramp.json"
    bc_filename =  folder * "bc_ramp.json"
    params_filename = folder * "params_ramp.json"

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

    params_data = parse_json(params_file)
    params_data["simulation_params"]["Initial time"] = t_initial
    params_data["simulation_params"]["Final time"] = t_final

    open(params_filename, "w") do f
        JSON.print(f, params_data, 2)
    end
end 
