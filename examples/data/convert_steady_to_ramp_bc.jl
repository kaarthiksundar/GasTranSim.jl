
using JSON
using GasTranSim


bc_ramp = Dict{String, Any}(
"boundary_pslack" => Dict(), 
"boundary_nonslack_flow" => Dict(), 
"boundary_compressor" => Dict())

t_initial = 0
t_ramp = 43200
t_final = 86400
p_init = 1e5

bc_steady_file = "./8-node-ramp-steady/bc_steady.json"
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
    bc_ramp["boundary_compressor"][key]["value"] = [val, val, val]
end


filename = "./8-node-ramp-steady/bc_ramp.json"
open(filename, "w") do f 
    JSON.print(f, bc_ramp, 2)
end


network_file = "./8-node-ramp-steady/network.json"
network_data  = parse_json(network_file)

ic_filename = "./8-node-ramp-steady/ic_ramp.json"

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



