using GasTranSim
using GLMakie

save_figures = true

base_path = split(Base.active_project(), "Project.toml")[1]
folder = base_path * "data/8-node/"
output_plot = base_path * "output/plots/"
tmp = base_path * "tmp/"

# eos = :ideal
eos = :simple_cnga
methods = [
    :explicit_staggered_grid,
    # :explicit_staggered_grid_new,
    :implicit_parabolic,
    # :explicit_hyperbolic,
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
    pipe_segs = nothing
    if method == :implicit_hyperbolic
        pipe_segs = 50
    end

    @info "Running method $method"
    ts = initialize_simulator(folder; method = method, pipe_segments=pipe_segs, eos = eos)

    if method == :implicit_parabolic
        ts.params[:base_dt] = 100 * ts.params[:base_dt]
    end
    if method == :implicit_hyperbolic
        ts.params[:base_dt] = 50 * ts.params[:base_dt]
    end

    run_simulator!(ts; method = method, showprogress = true)

    t = ts.sol["time_points"] ./ 3600.0
    pressure_node_5 = ts.sol["nodes"]["5"]["pressure"] ./ 1e6
    in_flow_node_6 = ts.sol["pipes"]["1"]["in_flow"]
    out_flow_node_5 = ts.sol["pipes"]["5"]["out_flow"]
    out_flow_node_3 = ts.sol["pipes"]["2"]["out_flow"] .- ts.sol["pipes"]["3"]["in_flow"]
    cratio_1 = ts.sol["compressors"]["1"]["compression_ratio"]
    cratio_2 = ts.sol["compressors"]["2"]["compression_ratio"]
    cratio_3 = ts.sol["compressors"]["3"]["compression_ratio"]

    results[method] = Dict(
        :time => t,
        :in_flow_node_6 => in_flow_node_6,
        :pressure_node_5 => pressure_node_5,
        :out_flow_node_5 => out_flow_node_5,
        :out_flow_node_3 => out_flow_node_3,
        :cratio_1 => cratio_1,
        :cratio_2 => cratio_2,
        :cratio_3 => cratio_3,
    )
end

function padded_range(vals...; pad_frac = 0.05)
    lo = minimum(vcat(vals...))
    hi = maximum(vcat(vals...))
    span = max(hi - lo, eps(Float64))
    pad = pad_frac * span
    return (lo - pad, hi + pad)
end

ylims_inflow = padded_range([results[m][:in_flow_node_6] for m in methods]...)
ylims_pressure = padded_range([results[m][:pressure_node_5] for m in methods]...)
ylims_out5 = padded_range([results[m][:out_flow_node_5] for m in methods]...)
ylims_out3 = padded_range([results[m][:out_flow_node_3] for m in methods]...)
ylims_cr1 = padded_range([results[m][:cratio_1] for m in methods]...)
ylims_cr2 = padded_range([results[m][:cratio_2] for m in methods]...)
ylims_cr3 = padded_range([results[m][:cratio_3] for m in methods]...)

function plot_panel!(ax, metric::Symbol, methods, results, method_colors, method_linestyles, method_labels)
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

update_theme!(fonts = (; regular = "Helvetica", bold = "Helvetica bold"))
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (1100, 1200))
g = f[1, 1] = GridLayout()

ax1 = Axis(
    g[1, 1],
    xticks = [0, 6, 12, 18, 24],
    title = rich("In-flow at pipe #1 (kgs", superscript("-1"), ")"),
    limits = (nothing, nothing, ylims_inflow[1], ylims_inflow[2]),
)
plot_panel!(ax1, :in_flow_node_6, methods, results, method_colors, method_linestyles, method_labels)

ax2 = Axis(
    g[2, 1],
    xticks = [0, 6, 12, 18, 24],
    title = "Pressure at node #5 (MPa)",
    limits = (nothing, nothing, ylims_pressure[1], ylims_pressure[2]),
)
plot_panel!(ax2, :pressure_node_5, methods, results, method_colors, method_linestyles, method_labels)

ax3 = Axis(
    g[3, 1],
    xticks = [0, 6, 12, 18, 24],
    title = rich("Withdrawal at node #5 (kgs", superscript("-1"), ")"),
    limits = (nothing, nothing, ylims_out5[1], ylims_out5[2]),
)
plot_panel!(ax3, :out_flow_node_5, methods, results, method_colors, method_linestyles, method_labels)

ax4 = Axis(
    g[4, 1],
    xticks = [0, 6, 12, 18, 24],
    title = rich("Withdrawal at node #3 (kgs", superscript("-1"), ")"),
    limits = (nothing, nothing, ylims_out3[1], ylims_out3[2]),
)
plot_panel!(ax4, :out_flow_node_3, methods, results, method_colors, method_linestyles, method_labels)

ax5 = Axis(
    g[5, 1],
    xticks = [0, 6, 12, 18, 24],
    title = "Compressor ratio #1",
    limits = (nothing, nothing, ylims_cr1[1], ylims_cr1[2]),
)
plot_panel!(ax5, :cratio_1, methods, results, method_colors, method_linestyles, method_labels)

ax6 = Axis(
    g[6, 1],
    xticks = [0, 6, 12, 18, 24],
    title = "Compressor ratio #2",
    limits = (nothing, nothing, ylims_cr2[1], ylims_cr2[2]),
)
plot_panel!(ax6, :cratio_2, methods, results, method_colors, method_linestyles, method_labels)

ax7 = Axis(
    g[7, 1],
    xticks = [0, 6, 12, 18, 24],
    title = "Compressor ratio #3",
    xlabel = "time (hrs.)",
    limits = (nothing, nothing, ylims_cr3[1], ylims_cr3[2]),
)
plot_panel!(ax7, :cratio_3, methods, results, method_colors, method_linestyles, method_labels)

resize_to_layout!(f)

save_figures && save(output_plot * "8-node-network-comparison-$eos.png", f)
save(tmp * "8-node-network-comparison-$eos.png", f)
