using GasTranSim
import PyPlot
const plt = PyPlot
using LaTeXStrings
using TerminalExtensions




folder = "test/data/model8/"

##=== Run from t=0 to t=t1 ===##

ts = initialize_simulator(folder)
run_simulator!(ts)
write_output(ts; output_path = folder, 
	output_file = "output_steady.json", final_state_file = "ic_delete.json")

sol = parse_json(folder*"output_steady.json")
t = sol["time_points"]/3600 #hrs


for i in keys(get(sol, "nodes", []))
	pr_node = sol["nodes"][i]["pressure"]
	error = maximum( abs.(pr_node .- pr_node[1]) )/pr_node[1]
	@show error
end

println("run complete")




