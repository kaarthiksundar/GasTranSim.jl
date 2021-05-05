using JSON
using Dierckx
using ProgressMeter
import PyPlot
# using TerminalExtensions


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

# Temp set to 239.11 K instead of 288 as given to get given sound speed
folder = "./data/model1pipe_fast_transients/"

ts = initialize_simulator(folder)

run_simulator!(ts)


inlet_pr = ts.sol["nodes"]["1"]["pressure"]
inlet_density = ts.nominal_values[:density] * get_density(ts, inlet_pr / ts.nominal_values[:pressure])
outlet_pr = ts.sol["nodes"]["2"]["pressure"]
outlet_density = ts.nominal_values[:density] * get_density(ts, outlet_pr/ ts.nominal_values[:pressure] )

inlet_massflux = ts.sol["pipes"]["1"]["in_flow"]/ts.ref[:pipe][1]["area"]
inlet_vel = inlet_massflux ./ inlet_density
outlet_massflux = ts.sol["pipes"]["1"]["out_flow"]/ts.ref[:pipe][1]["area"]
outlet_vel = outlet_massflux ./ outlet_density

t = ts.sol["time_points"]/60 #mins


fig, ax = PyPlot.subplots(4, 2, figsize=(8, 6), sharex=true)
ax[1, 1].plot(t, inlet_pr/1e6)
ax[2, 1].plot(t, inlet_density)
ax[3, 1].plot(t, inlet_massflux)
ax[4, 1].plot(t, inlet_vel)

ax[1, 2].plot(t, outlet_pr/1e6)
ax[2, 2].plot(t, outlet_density)
ax[3, 2].plot(t, outlet_massflux)
ax[4, 2].plot(t, outlet_vel)


ax[4, 1].set_xlabel("time (mins)")
ax[4, 2].set_xlabel("time (mins)")
ax[1, 1].set_ylabel("pressure (MPa)")
ax[2, 1].set_ylabel("density (kg/m^3)")
ax[3, 1].set_ylabel("mass flux (kg / m^2 s)")
ax[4, 1].set_ylabel("vel (m/s)")
fig.tight_layout()
# display(fig)
# PyPlot.close(fig)
show()

