using NGTransientSimulator
import PyPlot
const plt = PyPlot
using LaTeXStrings
using TerminalExtensions




folder = "test/data/model8_paper_VG_AZ/"

##=== Run from t=0 to t=t1 ===##

ts = initialize_simulator(folder; case_name="start", case_types=[:params])
run_simulator!(ts)
write_output(ts; output_path = folder, 
	output_file = "output_restart_1.json", final_state_file = "ic_restart_1.json")

sol = parse_json(folder*"output_restart_1.json")
t = sol["time_points"]/3600 #hrs
dt1 = sol["time_step"]	
pr_node5_a = sol["nodes"]["5"]["pressure"]
inflow_node6_a = sol["pipes"]["1"]["in_flow"]
outflow_node2_a = sol["pipes"]["1"]["out_flow"]
println("first half run complete")


##=== Run from t=t1 to t=tf ===##


ts = initialize_simulator(folder; case_name="restart_1", 
	case_types=[:params, :ic])
run_simulator!(ts)
write_output(ts; output_path = folder, output_file = "output_restart_2.json")
println("second half run complete")

sol = parse_json(folder*"output_restart_2.json")
t1 = sol["time_points"]/3600  #hrs
dt2 = sol["time_step"]	
pr_node5_b = sol["nodes"]["5"]["pressure"]
inflow_node6_b = sol["pipes"]["1"]["in_flow"]
outflow_node2_b = sol["pipes"]["1"]["out_flow"]


##=== Concatenate results of first two runs ===##
t_cat = vcat(t, t1[2:end])
pr_node5_cat = vcat(pr_node5_a, pr_node5_b[2:end])
inflow_node6_cat = vcat(inflow_node6_a, inflow_node6_b[2:end])
outflow_node2_cat = vcat(outflow_node2_a, outflow_node2_b[2:end])

##=== Run from t=0 to t=tf ===##
ts = initialize_simulator(folder; case_name="full", case_types=[:params])
run_simulator!(ts)
println("Third full run complete")

t3 = ts.sol["time_points"]/3600 #hrs
pr_node5 = ts.sol["nodes"]["5"]["pressure"]
inflow_node6 = ts.sol["pipes"]["1"]["in_flow"]
outflow_node2 = ts.sol["pipes"]["1"]["out_flow"]


##=== Compute norm of difference ===##
e1 = sum(abs.(pr_node5_cat - pr_node5) ./ abs.(pr_node5)) / length(pr_node5)
e2 = sum(abs.(inflow_node6_cat - inflow_node6) ./ abs.(inflow_node6)) / length(inflow_node6)
e3 = sum(abs.(outflow_node2_cat - outflow_node2) ./ abs.(outflow_node2)) / length(outflow_node2)

@show e1, e2, e3
# if maximum([e1, e2, e3]) < eps
# 	println("Passed unit test")
# else
# 	@error "Failed unit test"
# end


##=== Plot results ===##

fig, ax = plt.subplots(2, 1, figsize=(6, 6), sharex=true)



ax[1, 1].plot(t, inflow_node6_a, "r.", t1,  inflow_node6_b, "g.", t3, inflow_node6, "k--", markevery=3)
ax[1, 1].legend(["First half", "Second half", "Full"]) #string interpolation
#ax[1,1].set_ylim([250, 340])

ax[2, 1].plot(t, pr_node5_a/1e6, "r.", t1, pr_node5_b/1e6, "g.", t3, pr_node5/1e6, "k--", markevery=3)
ax[2, 1].legend(["First half", "Second half", "Full"])
#ax[2,1].set_ylim([3, 5.5])


ax[2, 1].set_xlabel("time (hrs)")
ax[1, 1].set_ylabel(L"mass inflow at node 6  ($\mathrm{kg}\mathrm{s}^{-1}$)")
ax[2, 1].set_ylabel("pressure at node 5 (MPa)")
fig.tight_layout()



