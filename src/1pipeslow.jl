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

folder = "./data/model1pipe_slow_transients/"

ts = initialize_simulator(folder; eos = :ideal)

run_simulator!(ts)


inlet_pr = ts.sol["nodes"]["1"]["pressure"]
inlet_density = ts.nominal_values[:density] * get_density(ts, inlet_pr / ts.nominal_values[:pressure])
outlet_pr = ts.sol["nodes"]["2"]["pressure"]
outlet_density = ts.nominal_values[:density] * get_density(ts, outlet_pr/ ts.nominal_values[:pressure] )

inlet_massflux = ts.sol["pipes"]["1"]["in_flow"]/ts.ref[:pipe][1]["area"]
inlet_vel = inlet_massflux ./ inlet_density
outlet_massflux = ts.sol["pipes"]["1"]["out_flow"]/ts.ref[:pipe][1]["area"]
outlet_vel = outlet_massflux ./ outlet_density

t = ts.sol["time_points"]/3600 #hrs

ts = initialize_simulator(folder; eos = :simple_cnga)

run_simulator!(ts)


inlet_pr_cnga = ts.sol["nodes"]["1"]["pressure"]
inlet_density_cnga = ts.nominal_values[:density] * get_density(ts, inlet_pr_cnga / ts.nominal_values[:pressure])
outlet_pr_cnga = ts.sol["nodes"]["2"]["pressure"]
outlet_density_cnga = ts.nominal_values[:density] * get_density(ts, outlet_pr_cnga/ ts.nominal_values[:pressure] )

inlet_massflux_cnga = ts.sol["pipes"]["1"]["in_flow"]/ts.ref[:pipe][1]["area"]
inlet_vel_cnga = inlet_massflux_cnga ./ inlet_density_cnga
outlet_massflux_cnga = ts.sol["pipes"]["1"]["out_flow"]/ts.ref[:pipe][1]["area"]
outlet_vel_cnga = outlet_massflux_cnga ./ outlet_density_cnga



fig, ax = PyPlot.subplots(4, 2, figsize=(12, 6), sharex=true)
ax[1, 1].plot(t, inlet_pr/1e6, t, inlet_pr_cnga/1e6)
ax[1, 1].legend(["ideal", "non-ideal"])

ax[2, 1].plot(t, inlet_density, t, inlet_density_cnga)

ax[3, 1].plot(t, inlet_massflux, t, inlet_massflux_cnga)
ax[4, 1].plot(t, inlet_vel, t, inlet_vel_cnga)

ax[1, 2].plot(t, outlet_pr/1e6, t, outlet_pr_cnga/1e6)
ax[2, 2].plot(t, outlet_density, t, outlet_density_cnga)
ax[3, 2].plot(t, outlet_massflux, t, outlet_massflux_cnga)
ax[4, 2].plot(t, outlet_vel, t, outlet_vel_cnga)

ax[3, 2].set_ylim(239, 241)
ax[4, 1].set_xlabel("time (hrs)")
ax[4, 2].set_xlabel("time (hrs)")
ax[1, 1].set_ylabel("pressure (MPa)")
ax[2, 1].set_ylabel("density (kg/m^3)")
ax[3, 1].set_ylabel("mass flux (kg / m^2 s)")
ax[4, 1].set_ylabel("vel (m/s)")
fig.tight_layout()
# display(fig)
# PyPlot.close(fig)
show()

