using GasTranSim

base_path = split(Base.active_project(), "Project.toml")[1]
folder = base_path * "data/model30/"
output_json = base_path * "output/solution/"

##=== Run from t=0 to t=t1 ===##
ts = initialize_simulator(folder)
run_simulator!(ts)
write_output(ts; output_path = output_json, 
	output_file = "model30_output_steady_delete.json", final_state_file = "model30_ic_delete.json")

sol = parse_json(output_json * "model30_output_steady_delete.json")
t = sol["time_points"] / 3600.0 # hrs

for i in keys(get(sol, "nodes", []))
	pr_node = sol["nodes"][i]["pressure"]
	error = maximum( abs.(pr_node .- pr_node[1]) )/pr_node[1]
	@show error
end

@info "run complete"




