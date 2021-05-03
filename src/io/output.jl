# struct OutputData
#     initial_time::Float64 
#     final_time::Float64 
#     time_points_pressure::Vector{Float64}
#     time_points_flow_rate::Vector{Float64}
#     junction::Dict{Int64,Any}
#     pipe::Dict{Int64,Any}
#     compressor::Dict{Int64,Any}
#     is_si_units::Bool
#     is_english_units::Bool 
#     is_per_unit::Bool
# end 

function initialize_output_data(ts::TransientSimulator)::OutputData 
    initial_time = ref(ts, :current_time) 
    final_time = params(ts, :tf)
    time_points_pressure = Float64[]
    time_points_flow_rate = Float64[]
    push!(time_points_pressure, initial_time)
    push!(time_points_flow_rate, initial_time)
    junction = Dict{Int64,Any}()
    pipe = Dict{Int64,Any}()
    compressor = Dict{Int64,Any}()
    # for (i, node) in get(ref(ts), :node, [])
    #     junction[i] = Dict(
    #         "pressure" => [],
    #         "flow" => []
    #     )


end 

struct OutputIntermediate
    time_pressure::Vector{Float64}
    time_flux::Vector{Float64}
    junction::Dict{Int64, Vector{Float64}}
    pipe::Dict{Int64, Tuple{ Vector{Float64}, Vector{Float64} } }
end 

struct OutputData
    initial_time::Float64
    final_time::Float64
    junction::Dict{Int64, Any}
    pipe::Dict{Int64, Dict{Any, Any}}
end 


function initialize_output_struc(ts::TransientSimulator)::OutputIntermediate

	out_int = OutputIntermediate([], [], Dict(), Dict())
	for (i, dummy) in ts.ref[:node]
        out_int.junction[i] = []
    end

    for (i, dummy) in ts.ref[:pipe]
        out_int.pipe[i] = ([], [])
    end
    return out_int
end


function update_output_struc!(ts::TransientSimulator, out_int::OutputIntermediate)

	push!(out_int.time_pressure, ts.ref[:current_time])
	push!(out_int.time_flux, ts.ref[:current_time] - ts.params[:dt]/2)

	for (i, dummy) in ts.ref[:node]
        push!(out_int.junction[i], ts.ref[:node][i]["pressure"])
    end

    for (i, dummy) in ts.ref[:pipe]
        push!(out_int.pipe[i][1], ts.ref[:pipe][i]["fr_mass_flux"])
        push!(out_int.pipe[i][2], ts.ref[:pipe][i]["to_mass_flux"])

    end
    return
end

function create_output(ts::TransientSimulator, out_int::OutputIntermediate)::OutputData

	out = OutputData(ts.params[:t_0], ts.ref[:current_time], Dict(), Dict())

	for (i, dummy) in ts.ref[:node]
        out.junction[i] = Spline1D(out_int.time_pressure, out_int.junction[i], k=1)

    end

    for (i, dummy) in ts.ref[:pipe]
    	out.pipe[i] = Dict("fr_mass_flux"=>Function, "to_mass_flux"=>Function, "density_profile" => Function,
    		"mass_flux_profile" => Function)
        out.pipe[i]["fr_mass_flux"] = Spline1D(out_int.time_flux, out_int.pipe[i][1], k=1)
        out.pipe[i]["to_mass_flux"] = Spline1D(out_int.time_flux, out_int.pipe[i][2], k=1) 

        n = ts.ref[:pipe][i]["num_discretization_points"]
        dx = ts.ref[:pipe][i]["dx"]
        L = ts.ref[:pipe][i]["length"]
        x_rho = LinRange(0, 1, n)
        x_phi = x_rho .+ dx/(2*L)
        rho = ts.ref[:pipe][i]["density_profile"][2:n-1] #len n
        phi = ts.ref[:pipe][i]["mass_flux_profile"][2:n] #len n+1

        out.pipe[i]["density_profile"] = Spline1D(x_rho[2:n-1], rho, k=1)
        out.pipe[i]["mass_flux_profile"] = Spline1D(x_phi[1:n-1], phi, k=1)


    end
    return out

end