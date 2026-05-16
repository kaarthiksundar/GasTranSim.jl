using GasTranSim
using GLMakie

save_figures = true

base_path = split(Base.active_project(), "Project.toml")[1]
folder = base_path * "data/1-pipe-slow-transients/"
output_plot = base_path * "output/plots/"
tmp = base_path * "tmp/"

methods = [
    :explicit_staggered_grid,
    :explicit_staggered_grid_new,
    :implicit_parabolic,
    :explicit_hyperbolic,
    :implicit_hyperbolic,
]

method_labels = Dict(
    :explicit_staggered_grid => "explicit_staggered_grid",
    :explicit_staggered_grid_new => "explicit_staggered_grid_new",
    :implicit_parabolic => "implicit_parabolic",
    :explicit_hyperbolic => "explicit_hyperbolic",
    :implicit_hyperbolic => "implicit_hyperbolic",
)

method_colors = Dict(
    :explicit_staggered_grid => :royalblue,
    :explicit_staggered_grid_new => :seagreen,
    :implicit_parabolic => :darkorange,
    :explicit_hyperbolic => :firebrick,
    :implicit_hyperbolic => :purple4,
)

method_linestyles = Dict(
    :explicit_staggered_grid => :solid,
    :explicit_staggered_grid_new => :dash,
    :implicit_parabolic => :dot,
    :explicit_hyperbolic => :dashdot,
    :implicit_hyperbolic => :dashdotdot,
)

results = Dict{Symbol,Dict{Symbol,Vector{Float64}}}()

for method in methods
    @info "Running method $method"
    ts = initialize_simulator(folder; method = method, eos = :ideal)
    run_simulator!(ts; method = method)

    nominal_density = nominal_values(ts, :density)
    nominal_pressure = nominal_values(ts, :pressure)
    solution = ts.sol
    pipe_solution = solution["pipes"]
    node_solution = solution["nodes"]

    inlet_pressure = node_solution["1"]["pressure"]
    outlet_pressure = node_solution["2"]["pressure"]

    inlet_density = nominal_density .* get_density.(Ref(ts), inlet_pressure ./ nominal_pressure)
    outlet_density = nominal_density .* get_density.(Ref(ts), outlet_pressure ./ nominal_pressure)

    area = ref(ts, :pipe, 1, "area")
    inlet_mass_flux = pipe_solution["1"]["in_flow"] ./ area
    outlet_mass_flux = pipe_solution["1"]["out_flow"] ./ area

    inlet_velocity = inlet_mass_flux ./ inlet_density
    outlet_velocity = outlet_mass_flux ./ outlet_density

    t = solution["time_points"] ./ 3600.0

    results[method] = Dict(
        :time => t,
        :inlet_pressure => inlet_pressure ./ 1e6,
        :outlet_pressure => outlet_pressure ./ 1e6,
        :inlet_density => inlet_density,
        :outlet_density => outlet_density,
        :inlet_mass_flux => inlet_mass_flux,
        :outlet_mass_flux => outlet_mass_flux,
        :inlet_velocity => inlet_velocity,
        :outlet_velocity => outlet_velocity,
    )
end

function padded_range(vals...; pad_frac = 0.05)
    lo = minimum(vcat(vals...))
    hi = maximum(vcat(vals...))
    span = max(hi - lo, eps(Float64))
    pad = pad_frac * span
    return (lo - pad, hi + pad)
end

ylims_pressure = padded_range([results[m][:inlet_pressure] for m in methods]..., [results[m][:outlet_pressure] for m in methods]...)
ylims_density = padded_range([results[m][:inlet_density] for m in methods]..., [results[m][:outlet_density] for m in methods]...)
ylims_mass_flux = padded_range([results[m][:inlet_mass_flux] for m in methods]..., [results[m][:outlet_mass_flux] for m in methods]...)
ylims_velocity = padded_range([results[m][:inlet_velocity] for m in methods]..., [results[m][:outlet_velocity] for m in methods]...)

update_theme!(fonts = (; regular = "Helvetica", bold = "Helvetica bold"))
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (1200, 800))
ga = f[1, 1] = GridLayout()
gb = f[2, 1] = GridLayout()
gc = f[3, 1] = GridLayout()
gd = f[4, 1] = GridLayout()

function plot_panel!(ax, metric::Symbol)
    handles = Any[]
    labels = String[]
    for method in methods
        h = lines!(
            ax,
            results[method][:time],
            results[method][metric],
            color = method_colors[method],
            linestyle = method_linestyles[method],
            linewidth = 2,
        )
        push!(handles, h)
        push!(labels, method_labels[method])
    end
    axislegend(ax, handles, labels, position = :rt)
    return
end

ax1 = Axis(ga[1, 1], title = "Inlet pressure (MPa)", limits = (nothing, nothing, ylims_pressure[1], ylims_pressure[2]))
plot_panel!(ax1, :inlet_pressure)

ax2 = Axis(ga[1, 2], title = "Outlet pressure (MPa)", limits = (nothing, nothing, ylims_pressure[1], ylims_pressure[2]))
plot_panel!(ax2, :outlet_pressure)

ax3 = Axis(gb[1, 1], title = rich("Inlet density (kgm", superscript("-3"), ")"), limits = (nothing, nothing, ylims_density[1], ylims_density[2]))
plot_panel!(ax3, :inlet_density)

ax4 = Axis(gb[1, 2], title = rich("Outlet density (kgm", superscript("-3"), ")"), limits = (nothing, nothing, ylims_density[1], ylims_density[2]))
plot_panel!(ax4, :outlet_density)

ax5 = Axis(
    gc[1, 1],
    title = rich("Inlet mass flux (kgm", superscript("-2"), "s", superscript("-1"), ")"),
    limits = (nothing, nothing, ylims_mass_flux[1], ylims_mass_flux[2]),
)
plot_panel!(ax5, :inlet_mass_flux)

ax6 = Axis(
    gc[1, 2],
    title = rich("Outlet mass flux (kgm", superscript("-2"), "s", superscript("-1"), ")"),
    limits = (nothing, nothing, ylims_mass_flux[1], ylims_mass_flux[2]),
)
plot_panel!(ax6, :outlet_mass_flux)

ax7 = Axis(
    gd[1, 1],
    title = rich("Inlet velocity (ms", superscript("-1"), ")"),
    xlabel = "time (hrs.)",
    limits = (nothing, nothing, ylims_velocity[1], ylims_velocity[2]),
)
plot_panel!(ax7, :inlet_velocity)

ax8 = Axis(
    gd[1, 2],
    title = rich("Outlet velocity (ms", superscript("-1"), ")"),
    xlabel = "time (hrs.)",
    limits = (nothing, nothing, ylims_velocity[1], ylims_velocity[2]),
)
plot_panel!(ax8, :outlet_velocity)

(save_figures) && save(output_plot * "method-comparison-ideal.png", f)
save(tmp * "method-comparison-ideal.png", f)
