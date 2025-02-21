
using JSON
using GasTranSim

bc_steady = Dict{String, Any}(
"boundary_pslack" => Dict(), 
"boundary_nonslack_flow" => Dict(), 
"boundary_compressor" => Dict())

bc_file = "./elpaso/bc.json"
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


filename = "./elpaso/bc_steady.json"
open(filename, "w") do f 
    JSON.print(f, bc_steady, 2)
end




