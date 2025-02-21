
using GasTranSim
using JSON

# Inputs
folder_in = "./elpaso/"
network_file =  folder_in * "network.json"
bc_steady_file = folder_in * "bc_steady.json"
params_file = folder_in * "params.json"

t_initial = 0
t_ramp = 3600*6
t_final = 3600*24
p_init = 1e5

# Outputs
folder_out =  folder_in
ic_filename = folder_out * "ic_ramp.json"
bc_filename =  folder_out * "bc_ramp.json"
params_filename = folder_out * "params_ramp.json"

bc_ramp = Dict{String, Any}(
"boundary_pslack" => Dict(), 
"boundary_nonslack_flow" => Dict(), 
"boundary_compressor" => Dict())

bc_steady_data = parse_json(bc_steady_file)

for (key, p_val) in bc_steady_data["boundary_pslack"]
    global  p_init = p_val
    break
end

for (key, p_slack) in bc_steady_data["boundary_pslack"]
    bc_ramp["boundary_pslack"][key] = Dict{String, Any}()
    bc_ramp["boundary_pslack"][key]["time"] = [t_initial, t_ramp, t_final] 
    bc_ramp["boundary_pslack"][key]["value"] = [p_init, p_slack, p_slack]
end

for (key, flow_val) in bc_steady_data["boundary_nonslack_flow"]
    bc_ramp["boundary_nonslack_flow"][key] = Dict{String, Any}()
    bc_ramp["boundary_nonslack_flow"][key]["time"] = [t_initial, t_ramp, t_final] 
    bc_ramp["boundary_nonslack_flow"][key]["value"] = [0, flow_val, flow_val]
end

for (key, compressor_dict) in bc_steady_data["boundary_compressor"]
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


for (node_id, -) in network_data["nodes"]
    ic_data["nodal_pressure"][node_id] = p_init
end

for (pipe_id, -) in network_data["pipes"]
    ic_data["pipe_flow"][pipe_id] = 0
end

for (comp_id, -) in network_data["compressors"]
    ic_data["compressor_flow"][comp_id] = 0
end

open(ic_filename, "w") do f
    JSON.print(f, ic_data, 2)
end

params_data = parse_json(params_file)
params_data["simulation_params"]["Initial time"] =  t_initial
params_data["simulation_params"]["Final time"] =  t_final


open(params_filename, "w") do f
    JSON.print(f, params_data, 2)
end