
using GasTranSim
using GLMakie

save_figures = false

base_path = split(Base.active_project(), "Project.toml")[1]
folder = base_path * "data/1-pipe-slow-transients/"
output_plot = base_path * "output/plots/"
output_json = base_path * "output/solution/"
tmp = base_path * "tmp/"

""" 
    Note: Temp for ideal gas case set to 239.11 K, 
    nonideal 288.7K to match pressure and density
"""

# Ideal EoS run

ts = initialize_simulator(folder; eos = :ideal)
run_simulator!(ts)

@info "ideal run completed"

nominal_density = nominal_values(ts, :density)
nominal_pressure = nominal_values(ts, :pressure)
solution = ts.sol
pipe_solution = solution["pipes"]
node_solution = solution["nodes"]

inlet_pressure = node_solution["1"]["pressure"]
inlet_density = nominal_density * get_density(ts, inlet_pressure / nominal_pressure)
outlet_pressure = node_solution["2"]["pressure"]
outlet_density = nominal_density * get_density(ts, outlet_pressure / nominal_pressure)

area = ref(ts, :pipe, 1, "area")
inlet_mass_flux = pipe_solution["1"]["in_flow"] / area
inlet_velocity = inlet_mass_flux ./ inlet_density
outlet_mass_flux = pipe_solution["1"]["out_flow"] / area
outlet_velocity = outlet_mass_flux ./ outlet_density

t = ts.sol["time_points"] / 3600.0 # hrs 

# CNGA EoS run

ts_cnga = initialize_simulator(
    folder;
    eos = :simple_cnga,
    case_name = "cnga",
    case_types = [:params],
)
run_simulator!(ts_cnga)
@info "simple CNGA run completed"

t_cnga = ts_cnga.sol["time_points"] / 3600.0 # hrs

nominal_density = nominal_values(ts_cnga, :density)
nominal_pressure = nominal_values(ts_cnga, :pressure)
solution = ts_cnga.sol
pipe_solution = solution["pipes"]
node_solution = solution["nodes"]

inlet_pressure_cnga = node_solution["1"]["pressure"]
inlet_density_cnga =
    nominal_density * get_density(ts_cnga, inlet_pressure_cnga / nominal_pressure)
outlet_pressure_cnga = node_solution["2"]["pressure"]
outlet_density_cnga =
    nominal_density * get_density(ts_cnga, outlet_pressure_cnga / nominal_pressure)

area = ref(ts, :pipe, 1, "area")
inlet_mass_flux_cnga = pipe_solution["1"]["in_flow"] / area
inlet_velocity_cnga = inlet_mass_flux_cnga ./ inlet_density_cnga
outlet_mass_flux_cnga = pipe_solution["1"]["out_flow"] / area
outlet_velocity_cnga = outlet_mass_flux_cnga ./ outlet_density_cnga

# Plot the results

function add_legend(ax, data, key)
    axislegend(ax, data, key, position = :rt, orientation = :horizontal)
end

update_theme!(fonts = (; regular = "Helvetica", bold = "Helvetica bold"))
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (1000, 700))
ga = f[1, 1] = GridLayout()
gb = f[2, 1] = GridLayout()
gc = f[3, 1] = GridLayout()
gd = f[4, 1] = GridLayout()

axmain = Axis(ga[1, 1], title = "Inlet pressure (MPa)")
ideal = scatterlines!(
    axmain,
    t,
    round.(inlet_pressure/1e6; digits = 2),
    alpha = 0.7,
    color = :orange,
    marker = :circle,
)
cnga = scatterlines!(
    axmain,
    t_cnga,
    round.(inlet_pressure_cnga/1e6; digits = 2),
    alpha = 0.7,
    linestyle = (:dashdot, :dense),
    color = :green,
    marker = :star4,
)
add_legend(axmain, [ideal, cnga], ["ideal", "non-ideal"])

axmain = Axis(ga[1, 2], title = "Outlet pressure (MPa)")
ideal = scatterlines!(
    axmain,
    t,
    round.(outlet_pressure/1e6; digits = 2),
    alpha = 0.7,
    color = :orange,
    marker = :circle,
)
cnga = scatterlines!(
    axmain,
    t_cnga,
    round.(outlet_pressure_cnga/1e6; digits = 2),
    alpha = 0.7,
    linestyle = (:dashdot, :dense),
    color = :green,
    marker = :star4,
)

axmain = Axis(gb[1, 1], title = rich("Inlet density (kgm", superscript("-3"), ")"))
ideal = scatterlines!(
    axmain,
    t,
    round.(inlet_density; digits = 2),
    alpha = 0.7,
    color = :orange,
    marker = :circle,
)
cnga = scatterlines!(
    axmain,
    t_cnga,
    round.(inlet_density_cnga; digits = 2),
    alpha = 0.7,
    linestyle = (:dashdot, :dense),
    color = :green,
    marker = :star4,
)

axmain = Axis(gb[1, 2], title = rich("Outlet density (kgm", superscript("-3"), ")"))
ideal = scatterlines!(
    axmain,
    t,
    round.(outlet_density; digits = 2),
    alpha = 0.7,
    color = :orange,
    marker = :circle,
)
cnga = scatterlines!(
    axmain,
    t_cnga,
    round.(outlet_density_cnga; digits = 2),
    alpha = 0.7,
    linestyle = (:dashdot, :dense),
    color = :green,
    marker = :star4,
)

axmain = Axis(
    gc[1, 1],
    title = rich("Inlet mass flux (kgm", superscript("-2"), "s", superscript("-1"), ")"),
    yticks = [150, 250, 350],
)
ideal = scatterlines!(
    axmain,
    t,
    round.(inlet_mass_flux; digits = 2),
    alpha = 0.7,
    color = :orange,
    marker = :circle,
)
cnga = scatterlines!(
    axmain,
    t_cnga,
    round.(inlet_mass_flux_cnga; digits = 2),
    alpha = 0.7,
    linestyle = (:dashdot, :dense),
    color = :green,
    marker = :star4,
)

axmain = Axis(
    gc[1, 2],
    title = rich("Outlet mass flux (kgm", superscript("-2"), "s", superscript("-1"), ")"),
    yticks = [100, 250, 400],
)
ideal = scatterlines!(
    axmain,
    t,
    round.(outlet_mass_flux; digits = 2),
    alpha = 0.7,
    color = :orange,
    marker = :circle,
)
cnga = scatterlines!(
    axmain,
    t_cnga,
    round.(outlet_mass_flux_cnga; digits = 2),
    alpha = 0.7,
    linestyle = (:dashdot, :dense),
    color = :green,
    marker = :star4,
)

axmain = Axis(
    gd[1, 1],
    title = rich("Inlet velocity (ms", superscript("-1"), ")"),
    xlabel = "time (hrs.)",
)
ideal = scatterlines!(
    axmain,
    t,
    round.(inlet_velocity; digits = 2),
    alpha = 0.7,
    color = :orange,
    marker = :circle,
)
cnga = scatterlines!(
    axmain,
    t_cnga,
    round.(inlet_velocity_cnga; digits = 2),
    alpha = 0.7,
    linestyle = (:dashdot, :dense),
    color = :green,
    marker = :star4,
)

axmain = Axis(
    gd[1, 2],
    title = rich("Outlet velocity (ms", superscript("-1"), ")"),
    xlabel = "time (hrs.)",
)
ideal = scatterlines!(
    axmain,
    t,
    round.(outlet_velocity; digits = 2),
    alpha = 0.7,
    color = :orange,
    marker = :circle,
)
cnga = scatterlines!(
    axmain,
    t_cnga,
    round.(outlet_velocity_cnga; digits = 2),
    alpha = 0.7,
    linestyle = (:dashdot, :dense),
    color = :green,
    marker = :star4,
)

(save_figures) && save(output_plot * "1-pipe-slow.png", f)
save(tmp * "1-pipe-slow.png", f)
