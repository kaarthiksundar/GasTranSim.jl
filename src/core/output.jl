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

function populate_solution!(ts::TransientSimulator, output::OutputData)
    dt = params(ts, :output_dt)
    dx = params(ts, :output_dx)
    times = collect(range(params(ts, :t_0), params(ts, :t_f), step=dt)) 
    data = ts.data
    sol = ts.sol
    units = params(ts, :units)
    units = 1

    function pressure_convertor(pu) 
        (units == 0) && (return pu * nominal_values(ts, :pressure)) 
        return pascal_to_psi(pu * nominal_values(ts, :pressure))
    end 

    function flow_convertor(pu)
        kgps_to_mmscfd = get_kgps_to_mmscfd_conversion_factor(params(ts))
        (units == 0) && (return pu * nominal_values(ts, :mass_flow)) 
        return pu * nominal_values(ts, :mass_flow) * kgps_to_mmscfd
    end 

    time_convertor(pu) = pu * nominal_values(ts, :time)
    
    function length_convertor(pu) 
        (units == 0) && (return pu * nominal_values(ts, :length))
        return m_to_miles(pu * nominal_values(ts, :length))
    end 

    
    for (i, _) in get(data, "nodes", [])
        key =  isa(i, String) ? parse(Int64, i) : i
        pressure_spl = output.node[key]["pressure"]
        pressure = [pressure_spl(t) for t in times]
        sol["nodes"][i]["pressure"] = pressure_convertor.(pressure)
        @show sol["nodes"][i]["pressure"]
    end

    for (i, _) in get(data, "pipes", [])
        key =  isa(i, String) ? parse(Int64, i) : i
        area = ref(ts, :pipe, key, "area")
        fr_flux_spl = output.pipe[key]["fr_mass_flux"]
        to_flux_spl = output.pipe[key]["to_mass_flux"]
        fr_flow = [fr_flux_spl(t) * area for t in times]
        to_flow = [to_flux_spl(t) * area for t in times]
        sol["pipes"][i]["in_flow"] = flow_convertor.(fr_flow)
        sol["pipes"][i]["out_flow"] = flow_convertor.(to_flow)
        fr_node = string(ref(ts, :pipe, key, "fr_node"))
        to_node = string(ref(ts, :pipe, key, "to_node"))
        sol["pipes"][i]["in_pressure"] = sol["nodes"][fr_node]["pressure"]
        sol["pipes"][i]["out_pressure"] = sol["nodes"][to_node]["pressure"]
    end 

    for (i, _) in get(data, "compressors", [])
        key =  isa(i, String) ? parse(Int64, i) : i
        fr_node = string(ref(ts, :compressor, key, "fr_node"))
        to_node = string(ref(ts, :compressor, key, "to_node"))
        sol["compressors"][i]["suction_pressure"] = sol["nodes"][fr_node]["pressure"]
        sol["compressors"][i]["discharge_pressure"] = sol["nodes"][to_node]["pressure"]
    end 

    sol["times"] = time_convertor.(times)

    for (i, value) in output.node 
        key =  isa(i, Int) ? string(i) : i
        pressure = value["pressure"]
    end 

    for (i, value) in output.pipe
        key =  isa(i, Int) ? string(i) : i
        L = ref(ts, :pipe, key, "length")
    end 
end 