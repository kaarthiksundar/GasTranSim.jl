mutable struct OutputData
    initial_time::Float64 
    final_time::Float64 
    node::Dict{Int64,Any}
    pipe::Dict{Int64,Any}
end 

function OutputData(ts::TransientSimulator)::OutputData 
    initial_time = ref(ts, :current_time)
    final_time = NaN 
    node = Dict{Int64,Any}()
    pipe = Dict{Int64,Any}()
    for (i, dummy) in ref(ts, :node)
        node[i] = Dict{String,Any}() 
    end 
    for (i, dummy) in ref(ts, :pipe)
        pipe[i] = Dict{String,Any}() 
    end 
    return OutputData(initial_time, final_time, node, pipe)
end 

struct OutputState
    time_pressure::Vector{Float64}
    time_flux::Vector{Float64}
    node::Dict{Int64,Any}
    pipe::Dict{Int64,Any}
end 

function initialize_output_state(ts::TransientSimulator)::OutputState 
    time_pressure = [ref(ts, :current_time)]
    time_flux = [ref(ts, :current_time)]
    node = Dict{Int64,Any}()
    pipe = Dict{Int64,Any}() 
    compressor = Dict{Int64,Any}()
    for i in keys(get(ref(ts), :node, []))
        node[i] = Dict(
            "pressure" => [ref(ts, :node, i, "initial_pressure")]
        )
    end 
    for i in keys(get(ref(ts), :pipe, []))
        pipe[i] = Dict(
            "fr_mass_flux" => [ref(ts, :pipe, i, "initial_mass_flux")],
            "to_mass_flux" => [ref(ts, :pipe, i, "initial_mass_flux")], 
        )
    end 
    return OutputState(time_pressure, time_flux, node, pipe)
end 

function update_output_state!(ts::TransientSimulator, state::OutputState)
    push!(state.time_pressure, ref(ts, :current_time))
    push!(state.time_flux, ref(ts, :current_time) - params(ts, :dt)/2)
	for (i, dummy) in ref(ts, :node)
        push!(state.node[i]["pressure"], ref(ts, :node, i, "pressure"))
    end
    for (i, dummy) in ref(ts, :pipe)
        push!(state.pipe[i]["fr_mass_flux"], ref(ts, :pipe, i, "fr_mass_flux"))
        push!(state.pipe[i]["to_mass_flux"], ref(ts, :pipe, i, "to_mass_flux"))
    end
    return
end 

function update_output_data!(ts::TransientSimulator, 
    state::OutputState, data::OutputData)
    data.initial_time = params(ts, :t_0)
    data.final_time = ref(ts, :current_time)
    for (i, dummy) in ref(ts, :node)
        data.node[i]["pressure"] = Spline1D(
            state.time_pressure, state.node[i]["pressure"], k=1
        )
    end 
    for (i, pipe) in ref(ts, :pipe)
        data.pipe[i]["fr_mass_flux"] = Spline1D(
            state.time_flux, state.pipe[i]["fr_mass_flux"], k=1
        )
        data.pipe[i]["to_mass_flux"] = Spline1D(
            state.time_flux, state.pipe[i]["to_mass_flux"], k=1
        )
        n = pipe["num_discretization_points"]
        dx = pipe["dx"]
        L = pipe["length"]
        x_rho = LinRange(0, L, n)
        x_mid = x_rho[1:n-1] .+ dx/2.0
        x_phi = [0, x_mid..., L]
        rho = pipe["density_profile"] # len n
        phi = pipe["mass_flux_profile"] # len n+1

        data.pipe[i]["final_density_profile"] = Spline1D(x_rho, rho, k=1)
        data.pipe[i]["final_mass_flux_profile"] = Spline1D(x_phi, phi, k=1)
    end 
end 