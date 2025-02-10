mutable struct OutputData
    initial_time::Float64 
    final_time::Float64 
    node::Dict{Int64,Any}
    pipe::Dict{Int64,Any}
    compressor::Dict{Int64,Any}
    final_state::Dict{Any,Any}
end 

function OutputData(ts::TransientSimulator)::OutputData 
    initial_time = ref(ts, :current_time)
    final_time = NaN 
    node = Dict{Int64,Any}()
    pipe = Dict{Int64,Any}()
    compressor = Dict{Int64,Any}()

    final_state = Dict{Any,Any}(
        "nodal_pressure" => Dict{Int64,Any}(),
        "pipe_flow" => Dict{Int64,Any}(),
        "pipe_pressure" => Dict{Int64,Any}(),
        "compressor_flow" => Dict{Int64, Any}()
    )
    for (i, _) in ref(ts, :node)
        node[i] = Dict{String,Any}() 
    end 
    for (i, _) in ref(ts, :pipe)
        pipe[i] = Dict{String,Any}() 
    end 
    # network need not contain compressors
    for (i, _) in get(ref(ts), :compressor, [])
        compressor[i] = Dict{String,Any}() 
    end 

    return OutputData(initial_time, final_time, node, pipe, compressor, final_state)
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

    if params(ts, :load_adjust) == true
        for i in keys(get(ref(ts), :node, []))
            node[i] = Dict(
                "pressure" => [ref(ts, :node, i, "pressure")], 
                "load_reduction" => [ref(ts, :node, i, "load_reduction")]
            )
        end
    else
        for i in keys(get(ref(ts), :node, []))
            node[i] = Dict(
                "pressure" => [ref(ts, :node, i, "pressure")]
            )
        end  
    end

    for i in keys(get(ref(ts), :pipe, []))
        mass_flux_profile = ref(ts, :pipe, i, "mass_flux_profile")
        pipe[i] = Dict(
            "fr_mass_flux" => [mass_flux_profile[1]],
            "to_mass_flux" => [mass_flux_profile[end]], 
        )
    end 
    for i in keys(get(ref(ts), :compressor, []))
        compressor[i] = Dict(
            "flow" => [ref(ts, :compressor, i, "flow")]
        )
    end 
    return OutputState(time_pressure, time_flux, node, pipe, compressor)
end 

function update_output_state!(ts::TransientSimulator, state::OutputState)
    push!(state.time_pressure, ref(ts, :current_time))
    push!(state.time_flux, ref(ts, :current_time) - params(ts, :dt)/2)
	for i in keys(get(ref(ts), :node, []))
        push!(state.node[i]["pressure"], ref(ts, :node, i, "pressure"))
    end

    if params(ts, :load_adjust) == true
        for i in keys(get(ref(ts), :node, []))
            push!(state.node[i]["load_reduction"], ref(ts, :node, i, "load_reduction"))
        end
    end

    for i in keys(get(ref(ts), :pipe, []))
        push!(state.pipe[i]["fr_mass_flux"], ref(ts, :pipe, i, "fr_mass_flux"))
        push!(state.pipe[i]["to_mass_flux"], ref(ts, :pipe, i, "to_mass_flux"))
    end
    for i in keys(get(ref(ts), :compressor, []))
        push!(state.compressor[i]["flow"], ref(ts, :compressor, i, "flow"))
    end
    return
end 

function update_output_data!(ts::TransientSimulator, 
    state::OutputState, data::OutputData)
    data.initial_time = params(ts, :t_0)
    data.final_time = ref(ts, :current_time)
    for (i, _) in ref(ts, :node)
        data.node[i]["pressure"] = Spline1D(
            state.time_pressure, state.node[i]["pressure"], k=1
        )
        data.final_state["nodal_pressure"][i] = ref(ts, :node, i, "pressure")
    end 

    if params(ts, :load_adjust) == true
        for (i, _) in ref(ts, :node)
            data.node[i]["load_reduction"] = Spline1D(
                state.time_flux, state.node[i]["load_reduction"], k=1
            )
        end 
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
        area = pipe["area"]
        x_rho = LinRange(0, L, n)
        x_mid = x_rho[1:n-1] .+ dx/2.0
        # this is what needs to be done to replicate
        # x_phi = [-(dx/2), x_mid..., L+(dx/2)] 
        x_phi = [0.0, x_mid..., L]
        rho = pipe["density_profile"] # len n
        pressure = [get_pressure(ts, val) for val in rho]
        phi = pipe["mass_flux_profile"] # len n+1
        flow = phi .* area

        data.final_state["pipe_flow"][i] = Spline1D(x_phi, flow, k=1)
        data.final_state["pipe_pressure"][i] = Spline1D(x_rho, pressure, k=1)
    end 
    for i in keys(get(ref(ts), :compressor, []))
        data.compressor[i]["flow"] = Spline1D(
            state.time_flux, state.compressor[i]["flow"], k=1
        )
        data.final_state["compressor_flow"][i] = ref(ts, :compressor, i, "flow")
    end 
end 

function populate_solution!(ts::TransientSimulator, output::OutputData)
    dt = params(ts, :output_dt)
    times = collect(range(params(ts, :t_0), params(ts, :t_f), step=dt)) 
    data = ts.data
    sol = ts.sol
    units = params(ts, :units)

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

    sol["final_state"]["time"] = time_convertor(ref(ts, :current_time))
    sol["final_state"]["initial_nodal_pressure"] = Dict{String,Any}()
    
    for (i, _) in get(data, "nodes", [])
        key =  isa(i, String) ? parse(Int64, i) : i
        pressure_spl = output.node[key]["pressure"]
        pressure = [pressure_spl(t) for t in times]
        sol["nodes"][i]["pressure"] = pressure_convertor.(pressure)
        sol["final_state"]["initial_nodal_pressure"][i] = 
            pressure_convertor(output.final_state["nodal_pressure"][key])
    end



    if params(ts, :load_adjust) == true
        for (i, _) in ref(ts, :node)
            sol["nodes"][string(i)]["load_reduction"] = Vector{Float64}()
            load_reduction_spl = output.node[i]["load_reduction"]
            load_reduction = [load_reduction_spl(t) for t in times]
            sol["nodes"][string(i)]["load_reduction"] = flow_convertor.(load_reduction)
        end
        sol["load_reduction_nodes"] = Vector{String}()
        for node_id in ref(ts, :load_reduction_nodes)
            push!(sol["load_reduction_nodes"], string(node_id))
        end
    end


    sol["final_state"]["initial_pipe_flow"] = Dict{String,Any}()
    sol["final_state"]["initial_pipe_pressure"] = Dict{String,Any}()

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
        flow_spl = output.final_state["pipe_flow"][key]
        pressure_spl = output.final_state["pipe_pressure"][key]
        sol["final_state"]["initial_pipe_flow"][i] = Dict{String,Any}(
            "distance" => length_convertor.(get_knots(flow_spl)), 
            "value" => flow_convertor.(get_coeffs(flow_spl))
        )
        sol["final_state"]["initial_pipe_pressure"][i] = Dict{String,Any}(
            "distance" => length_convertor.(get_knots(pressure_spl)), 
            "value" => pressure_convertor.(get_coeffs(pressure_spl))
        )
    end 
    sol["final_state"]["initial_compressor_flow"] = Dict{String,Any}()
    for (i, _) in get(data, "compressors", [])
        key =  isa(i, String) ? parse(Int64, i) : i
        fr_node = string(ref(ts, :compressor, key, "fr_node"))
        to_node = string(ref(ts, :compressor, key, "to_node"))
        sol["compressors"][i]["suction_pressure"] = sol["nodes"][fr_node]["pressure"]
        sol["compressors"][i]["discharge_pressure"] = sol["nodes"][to_node]["pressure"]
        sol["compressors"][i]["compression_ratio"] = sol["nodes"][to_node]["pressure"] ./ sol["nodes"][fr_node]["pressure"]
        flow_spl = output.compressor[key]["flow"]
        flow = [flow_spl(t) for t in times]
        sol["compressors"][i]["flow"] = flow_convertor.(flow)
        sol["final_state"]["initial_compressor_flow"][i] = 
            flow_convertor(output.final_state["compressor_flow"][key])
    end 

    sol["time_step"] = time_convertor(params(ts, :dt))
    sol["time_points"] = time_convertor.(times)
    return
end 