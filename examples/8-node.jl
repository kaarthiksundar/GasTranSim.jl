using GasTranSim
using CairoMakie

save_figures = false
save_output = false

base_path = split(Base.active_project(), "Project.toml")[1]
folder = base_path * "data/8-node/"
output_plot = base_path * "output/plots/"
output_json = base_path * "output/solution/"
tmp = base_path * "tmp/"

##=== Run case ===##
ts = initialize_simulator(folder)
run_simulator!(ts)
@info "First run completed"

t = ts.sol["time_points"] / 3600.0 # hrs
dt1 = params(ts, :dt) * nominal_values(ts, :time)

cratio_1 = ts.sol["compressors"]["1"]["compression_ratio"]
cratio_2 = ts.sol["compressors"]["2"]["compression_ratio"]
cratio_3 = ts.sol["compressors"]["3"]["compression_ratio"]

pressure_node_5 = ts.sol["nodes"]["5"]["pressure"]

in_flow_node_6 = ts.sol["pipes"]["1"]["in_flow"]
out_flow_node_2 = ts.sol["pipes"]["1"]["out_flow"]

out_flow_node_5 = ts.sol["pipes"]["5"]["out_flow"]
out_flow_node_3 =  ts.sol["pipes"]["2"]["out_flow"] .- ts.sol["pipes"]["3"]["in_flow"]

##=== Run  case ===##
ts = initialize_simulator(folder; eos = :ideal, case_name="fine_time_step", case_types=[:params])
run_simulator!(ts)
@info "Second run complete"

t1 = ts.sol["time_points"] / 3600.0 # hrs
dt2 = params(ts, :dt) * nominal_values(ts, :time)

pressure_node_5_fine = ts.sol["nodes"]["5"]["pressure"]

in_flow_node_6_fine = ts.sol["pipes"]["1"]["in_flow"]
out_flow_node_2_fine = ts.sol["pipes"]["1"]["out_flow"]

##=== Plot results ===##

with_theme(theme_latexfonts()) do
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
        size = (700, 700))
    ga = f[1, 1] = GridLayout()
    gb = f[2, 1] = GridLayout()
    gc = f[3, 1] = GridLayout()
    gd = f[4, 1] = GridLayout()

    ax1 = Axis(ga[1, 1], xticks = [0, 6, 12, 18, 24], 
        yticks = [260, 300, 340], title = 
        rich("In-flow at pipe #1 (kgs", superscript("-1"), ")"))
    ax2 = Axis(gb[1, 1], xticks = [0, 6, 12, 18, 24], title = "Pressure at node #5 (MPa)")
    ax3 = Axis(gc[1, 1], xticks = [0, 6, 12, 18, 24], 
        yticks = [120, 150, 180], title = 
        rich("Withdrawal at node (kgs", superscript("-1"), ")"))
    ax4 = Axis(gd[1, 1], xticks = [0, 6, 12, 18, 24], 
        title = "Compressor ratio", xlabel="time (hrs.)")
    resize_to_layout!(f)

    coarse = lines!(ax1, t, in_flow_node_6, alpha = 0.7, color=:orange)
    fine = lines!(ax1, t1, in_flow_node_6_fine, alpha = 0.7, linestyle=(:dashdot, :dense), color=:green)
    axislegend(ax1, [coarse, fine], [L"\Delta t = %$(dt1) ~\text{sec.}", L"\Delta t = %$(dt2) ~\text{sec.}"], position = :cb, 
        orientation = :horizontal)


    coarse = lines!(ax2, t, round.(pressure_node_5/1e6; digits=2), alpha = 0.7, color=:orange)
    fine = lines!(ax2, t, round.(pressure_node_5_fine/1e6; digits=2), alpha = 0.7, linestyle=(:dashdot, :dense), color=:green)
    axislegend(ax2, [coarse, fine], [L"\Delta t = %$(dt1) ~\text{sec.}", L"\Delta t = %$(dt2) ~\text{sec.}"], position = :ct, 
        orientation = :horizontal)

    withdrawal_5 = lines!(ax3, t, out_flow_node_5, alpha = 0.7, color=:orange)
    withdrawal_3 = lines!(ax3, t, out_flow_node_3, alpha = 0.7, linestyle=(:dashdot, :dense), color=:green)
    axislegend(ax3, [withdrawal_5, withdrawal_3], ["#5", "#3"], position = :rt, orientation = :horizontal)

    r_1 = lines!(ax4, t, cratio_1, alpha = 0.7, color=:orange)
    r_2 = lines!(ax4, t, cratio_2, alpha = 0.7, linestyle=(:dashdot, :dense), color=:green)
    r_3 = lines!(ax4, t, cratio_3, alpha = 0.7, linestyle=(:dash, :dense), color=:blue)
    axislegend(ax4, [r_1, r_2, r_3], ["#1", "#2", "#3"], position = :rt, orientation = :horizontal)


    save_figures && save(output_plot * "8-node.png", f)
    save(tmp * "8-node.png", f)
end 

