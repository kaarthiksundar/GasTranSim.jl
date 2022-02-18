using GasTranSim
import PyPlot
const plt = PyPlot
using LaTeXStrings
using TerminalExtensions


folder = "examples/data/8-node/"

##=== Run case ===##

ts = initialize_simulator(folder)

run_simulator!(ts)

println("first run complete")


t = ts.sol["time_points"]/3600 #hrs
dt1 = ts.params[:dt] * ts.nominal_values[:time]	

cmp1 = ts.sol["compressors"]["1"]["compression_ratio"]
cmp2 = ts.sol["compressors"]["2"]["compression_ratio"]
cmp3 = ts.sol["compressors"]["3"]["compression_ratio"]

pr_node5 = ts.sol["nodes"]["5"]["pressure"]

inflow_node6 = ts.sol["pipes"]["1"]["in_flow"]
outflow_node2 = ts.sol["pipes"]["1"]["out_flow"]

outflow_node5 = ts.sol["pipes"]["5"]["out_flow"]
outflow_node3 =  ts.sol["pipes"]["2"]["out_flow"] .- ts.sol["pipes"]["3"]["in_flow"]


##=== Run  case ===##


ts = initialize_simulator(folder; eos = :ideal, case_name="fine_time_step", case_types=[:params])

run_simulator!(ts)

println("second run complete")

t1 = ts.sol["time_points"]/3600 #hrs
dt2 = ts.params[:dt] * ts.nominal_values[:time]	


pr_node5_fine = ts.sol["nodes"]["5"]["pressure"]

inflow_node6_fine = ts.sol["pipes"]["1"]["in_flow"]
outflow_node2_fine = ts.sol["pipes"]["1"]["out_flow"]


##=== Plot results ===##

fig, ax = plt.subplots(4, 1, figsize=(6, 12), sharex=true)



ax[1, 1].plot(t, inflow_node6)
ax[1,1].plot(t1, inflow_node6_fine)
ax[1, 1].legend(["pipe 1 inflow dt=$dt1", "pipe 1 inflow  dt=$dt2"]) #string interpolation
#ax[1,1].set_ylim([250, 340])

ax[2, 1].plot(t, pr_node5/1e6, t1, pr_node5_fine/1e6)
ax[2, 1].legend(["dt = $dt1", "dt=$dt2"])
#ax[2,1].set_ylim([3, 5.5])

ax[3, 1].plot(t, outflow_node5)
ax[3, 1].plot(t, outflow_node3)
ax[3, 1].legend(["node 5", "node 3"])

ax[4, 1].plot(t, cmp1)
ax[4, 1].plot(t, cmp2)
ax[4, 1].plot(t, cmp3)
ax[4, 1].legend(["comp 1", "comp 2", "comp 3"])


ax[4, 1].set_xlabel("time (hrs)")
ax[1, 1].set_ylabel(L"mass flow  ($\mathrm{kg}\mathrm{s}^{-1}$)")
ax[2, 1].set_ylabel("pressure at node 5 (MPa)")
ax[3, 1].set_ylabel(L"withdrawal ($\mathrm{kg}\mathrm{s}^{-1}$)")
ax[4, 1].set_ylabel("Compressor Ratios")
fig.tight_layout()


