using JSON
using Dierckx
using ProgressMeter
import PyPlot
using TerminalExtensions


include("io/json.jl")
include("io/data_utils.jl")

include("unit_conversion/unit_convertor_utils.jl")
include("unit_conversion/to_si.jl")
include("unit_conversion/to_english.jl")
include("unit_conversion/to_pu.jl")
include("unit_conversion/unit_convertors.jl")

include("core/eos.jl")
include("core/types.jl")
include("core/ref.jl")
include("core/bc.jl")
include("core/sol.jl")
include("core/initialize_ts.jl")
include("core/run_task.jl")
include("core/time_integration.jl")
include("core/run_ts.jl")
include("core/output.jl")

folder = "./data/model8_paper_VG_AZ/"

ts = initialize_simulator(folder)

run_simulator!(ts)

cmp1 = ts.sol["compressors"]["1"]["compression_ratio"]
cmp2 = ts.sol["compressors"]["2"]["compression_ratio"]
cmp3 = ts.sol["compressors"]["3"]["compression_ratio"]

pr_node5 = ts.sol["nodes"]["5"]["pressure"]

inflow_node1 = ts.sol["pipes"]["1"]["in_flow"]

outflow_node5 = ts.sol["pipes"]["5"]["out_flow"]
outflow_node3 =  ts.sol["pipes"]["2"]["out_flow"] .- ts.sol["pipes"]["3"]["in_flow"]

t = ts.sol["time_points"]/3600 #hrs


fig, ax = PyPlot.subplots(4, 1, figsize=(6, 12), sharex=true)
ax[2, 1].plot(t, pr_node5/1e6)
ax[1, 1].plot(t, inflow_node1)
ax[3, 1].plot(t, outflow_node5)
ax[3, 1].plot(t, outflow_node3)
ax[4, 1].plot(t, cmp1)
ax[4, 1].plot(t, cmp2)
ax[4, 1].plot(t, cmp3)
ax[4, 1].legend(["comp 1", "comp 2", "comp 3"])


ax[2, 1].set_xlabel("time (hrs)")

ax[2, 1].set_ylabel("pressure at node 5 (MPa)")
ax[1, 1].set_ylabel("mass inflow at node 1 (kg/s)")
ax[3, 1].set_ylabel("withdrawal at node 5 and 3 (kg/s)")
ax[4, 1].set_ylabel("Compressor Ratios")

fig.tight_layout()
# display(fig)
# PyPlot.close(fig)
show()

