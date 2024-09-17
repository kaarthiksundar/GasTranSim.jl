
using GasTranSim
import PyPlot
const plt = PyPlot
using LaTeXStrings
using TerminalExtensions




# Temp for ideal gas case set to 239.11 K, nonideal 288.7K to match p, rho
folder = "./data/model1pipe_fast_transients/"

##=== Run ideal gas case ===##

ts = initialize_simulator(folder; eos = :ideal)

run_simulator!(ts)

println("ideal completed")


inlet_pr = ts.sol["nodes"]["1"]["pressure"]
inlet_density = ts.nominal_values[:density] * get_density(ts, inlet_pr / ts.nominal_values[:pressure])
outlet_pr = ts.sol["nodes"]["2"]["pressure"]
outlet_density = ts.nominal_values[:density] * get_density(ts, outlet_pr/ ts.nominal_values[:pressure] )

inlet_massflux = ts.sol["pipes"]["1"]["in_flow"]/ts.ref[:pipe][1]["area"]
inlet_vel = inlet_massflux ./ inlet_density
outlet_massflux = ts.sol["pipes"]["1"]["out_flow"]/ts.ref[:pipe][1]["area"]
outlet_vel = outlet_massflux ./ outlet_density

t = ts.sol["time_points"]/60 #mins

##=== Run non-ideal gas case ===##

ts1 = initialize_simulator(folder; eos = :simple_cnga, case_name="cnga", case_types=[:params])

run_simulator!(ts1)

println("simple cnga completed")

t1 = ts1.sol["time_points"]/60 #mins

inlet_pr_cnga = ts1.sol["nodes"]["1"]["pressure"]
inlet_density_cnga = ts1.nominal_values[:density] * get_density(ts1, inlet_pr_cnga / ts1.nominal_values[:pressure])
outlet_pr_cnga = ts1.sol["nodes"]["2"]["pressure"]
outlet_density_cnga = ts1.nominal_values[:density] * get_density(ts1, outlet_pr_cnga/ ts1.nominal_values[:pressure] )

inlet_massflux_cnga = ts1.sol["pipes"]["1"]["in_flow"]/ts1.ref[:pipe][1]["area"]
inlet_vel_cnga = inlet_massflux_cnga ./ inlet_density_cnga
outlet_massflux_cnga = ts1.sol["pipes"]["1"]["out_flow"]/ts1.ref[:pipe][1]["area"]
outlet_vel_cnga = outlet_massflux_cnga ./ outlet_density_cnga


##=== Plot results ===##

fig, ax = plt.subplots(4, 2, figsize=(12, 6), sharex=true)
ax[1, 1].plot(t, inlet_pr/1e6, t1, inlet_pr_cnga/1e6)
ax[1, 1].legend(["ideal", "non-ideal"])

ax[2, 1].plot(t, inlet_density, t1, inlet_density_cnga)
ax[2, 1].set_ylim(55, 58)


ax[3, 1].plot(t, inlet_massflux, t1, inlet_massflux_cnga)
ax[4, 1].plot(t, inlet_vel, t1, inlet_vel_cnga)

ax[1, 2].plot(t, outlet_pr/1e6, t1, outlet_pr_cnga/1e6)
ax[2, 2].plot(t, outlet_density, t1, outlet_density_cnga)
ax[3, 2].plot(t, outlet_massflux, t1, outlet_massflux_cnga)
ax[4, 2].plot(t, outlet_vel, t1, outlet_vel_cnga)

ax[4, 1].set_xlabel("time (hrs)")
ax[4, 2].set_xlabel("time (hrs)")
ax[1, 1].set_ylabel("pressure (MPa)")
ax[2, 1].set_ylabel(L"$\rho \;(\mathrm{kg}\mathrm{m}^{-3})$")
ax[3, 1].set_ylabel(L"$\phi \; (\mathrm{kg}\mathrm{m}^{-2} \mathrm{s})$")
ax[4, 1].set_ylabel(L"$v \; (\mathrm{m}\mathrm{s}^{-1})$")
fig.tight_layout()

