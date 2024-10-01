using GLMakie
using JSON 

pascal_to_psi(pascal) = pascal / 6894.75729

pressure = Dict()
flow = Dict()

base_path = split(Base.active_project(), "Project.toml")[1]
output_anim = base_path * "output/animation/"
output_json = base_path * "output/solution/"
tmp = base_path * "tmp/"

first_file = output_json * "1-pipe-fast-state-1-min.json"
first_data = JSON.parsefile(first_file)
px = first_data["initial_pipe_pressure"]["1"]["distance"]
fx = first_data["initial_pipe_flow"]["1"]["distance"]

for i in 1:60
    file = output_json * "1-pipe-fast-state-$i-min.json"
    data = JSON.parsefile(file)
    pressure[i] = data["initial_pipe_pressure"]["1"]["value"]
    flow[i] = data["initial_pipe_flow"]["1"]["value"]
end

time = Observable(1)
p = @lift(round.(pascal_to_psi(pressure[$time]); digits=2))
f = @lift(round.(flow[$time]; digits=2))
f_control = Observable(Point2f[(20.0, flow[1][end])])
update_theme!(fonts = (; regular = "Helvetica", bold = "Helvetica bold"))
fig = Figure(fontsize = 12)
ax1 = Axis(fig[1, 1], yticklabelcolor = :blue, 
    ylabelcolor = :blue,
    title = @lift("Pressure and flow profiles in the pipe at t = $(round($time, digits = 1)) min.\n Pressure (flow) is controlled at the start (end) of pipe (denoted by dots)"), 
    xlabel = "distance from start of pipe (km)", 
    ylabel = "Pressure (psi)")
ax2 = Axis(fig[1, 1], yticklabelcolor = :red, yaxisposition = :right,
    ylabelcolor = :red, 
    ylabel = rich("Flow (kgs", superscript("-1"), ")")) 
slack_pressure = Observable(Point2f[(0.0, pascal_to_psi(
    first_data["initial_pipe_pressure"]["1"]["value"][1])
    )])
scatter!(ax1, slack_pressure, color = :blue)
scatter!(ax2, f_control, color = :red)
pressure_plot = lines!(ax1, round.(px/1000; digits=2), p, color = :blue, linewidth = 1.5)
flow_plot = lines!(ax2, round.(fx/1000.0; digits=2), f, color = :red, linewidth = 1.5)
ylims!(ax1, 400, 1000)
ylims!(ax2, -50, 850)

framerate = 3
timestamps = range(1, 60)

record(fig, output_anim * "1-pipe-fast-pipe-animation.mp4", timestamps;
        framerate = framerate) do t
    time[] = t
    f_control[] = Point2f[(20.0, flow[t][end])]
end