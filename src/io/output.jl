struct OutputData
    initial_time::Float64 
    final_time::Float64 
    node::Dict{Int64,Any}
    pipe::Dict{Int64,Any}
    compressor::Dict{Int64,Any}
end 

function initialize_output_data(ts::TransientSimulator)::OutputData 
    initial_time = ref(ts, :current_time) 
    final_time = ref(ts, :current_time) 
    times = [initial_time]
    node = Dict{Int64,Any}()
    pipe = Dict{Int64,Any}()
    compressor = Dict{Int64,Any}()
    for i in keys(get(ref(ts), :node, []))
        pressure = [ref(ts, :node, i, "initial_pressure")]
        withdrawal = [(-1.0) * ref(ts, :node, i, "initial_injection")]
        node[i] = Dict(
            "pressure" => Spline1D(times, pressure, k=1), 
            "withdrawal" => Spline1D(times, withdrawal, k=1)
        )
    end 
    for i in keys(get(ref(ts), :pipe, []))
        mass_flux = [ref(ts, :pipe, i, "initial_mass_flux")]
        n = ref(ts, :pipe, i, "num_discretization_points")
        dx = ref(ts, :pipe, i, "dx")
        L = ref(ts, :pipe, i, "length")
        x_density = LinRange(0, 1, n)
        x_mass_flux = x_density .+ dx/(2*L)
        density = ts.ref[:pipe][i]["density_profile"][2:n-1] #len n
        mass_flux = ts.ref[:pipe][i]["mass_flux_profile"][2:n] #len n+1
        pipe[i] = Dict(
            "final_density_profile" => Spline1D(x_density[2:n-1], density, k=1),
            "fr_mass_flux" => Spline(times, mass_flux, k=1), 
            "to_mass_flux" => Spline(times, mass_flux, k=1), 
            "final_mass_flux_profile" => Spline1D(x_mass_flux[2:n-1], mass_flux, k=1),
        )
    end 

    for i in keys(get(ref(ts), :compressor, []))
        c_ratio = [ref(ts, :compressor, i, "initial_c_ratio")]
        compressor[i] = Dict(
            "c_ratio" => Spline1D(times, c_ratio, k=1)
        )
    end 
    return OutputData(initial_time, final_time, node, pipe, compressor)
end 

struct OutputState
    time_pressure::Vector{Float64}
    time_flux::Vector{Float64}
    node::Dict{Int64,Any}
    pipe::Dict{Int64,Any}
    compressor::Dict{Int64,Any} 
end 


function initialize_output_state(ts::TransientSimulator)::OutputState 

    time_pressure = [ref(ts, :current_time)]
    time_flux = [ref(ts, :current_time)]
    node = Dict{Int64,Any}()
    pipe = Dict{Int64,Any}() 
    compressor = Dict{Int64,Any}()
    for i in keys(get(ref(ts), :node, []))
        pressure = [ref(ts, :node, i, "initial_pressure")]
        withdrawal = [(-1.0) * ref(ts, :node, i, "initial_injection")]
        node[i] = Dict(
            "pressure" => [pressure], 
            "withdrawal" => [withdrawal]
        )
    end 
    for i in keys(get(ref(ts), :pipe, []))
        mass_flux = [ref(ts, :pipe, i, "initial_mass_flux")]
        n = ref(ts, :pipe, i, "num_discretization_points")
        dx = ref(ts, :pipe, i, "dx")
        L = ref(ts, :pipe, i, "length")
        x_density = LinRange(0, 1, n)
        x_mass_flux = x_density .+ dx/(2*L)
        density = ts.ref[:pipe][i]["density_profile"][2:n-1] #len n
        mass_flux = ts.ref[:pipe][i]["mass_flux_profile"][2:n] #len n+1
        pipe[i] = Dict(
            "internal_density_profile" => Spline1D(x_density[2:n-1], density, k=1),
            "fr_mass_flux" => [mass_flux],
            "to_mass_flux" => [mass_flux], 
            "internal_mass_flux_profile" => Spline1D(x_mass_flux[1:n-1], mass_flux, k=1),
        )
    end 

    for i in keys(get(ref(ts), :compressor, []))
        c_ratio = [ref(ts, :compressor, i, "initial_c_ratio")]
        compressor[i] = Dict(
            "c_ratio" => [c_ratio]
        )
    end 

    return OutputState(time_pressure, time_flux, node, pipe, compressor)

end 

struct OutputIntermediate
    time_pressure::Vector{Float64}
    time_flux::Vector{Float64}
    junction::Dict{Int64, Vector{Float64}}
    pipe::Dict{Int64, Tuple{ Vector{Float64}, Vector{Float64} } }
end 

# struct OutputData
#     initial_time::Float64
#     final_time::Float64
#     junction::Dict{Int64, Any}
#     pipe::Dict{Int64, Dict{Any, Any}}
# end 


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