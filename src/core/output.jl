
function OutputData(ts::TransientSimulator)::OutputData
    initial_time = ref(ts, :current_time)
    final_time = NaN
    time_points = Vector{Float64}()
    node = Dict{Int64,Any}()
    pipe = Dict{Int64,Any}()
    compressor = Dict{Int64,Any}()

    final_state = Dict{Any,Any}(
        "nodal_pressure" => Dict{Int64,Any}(),
        "pipe_flow" => Dict{Int64,Any}(),
        "pipe_pressure" => Dict{Int64,Any}(),
        "compressor_flow" => Dict{Int64,Any}(),
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

    return OutputData(initial_time, final_time, time_points, node, pipe, compressor, final_state)
end

function initialize_output_state(ts::TransientSimulator)::OutputState
    current_time = ref(ts, :current_time)
    output_dt = params(ts, :output_dt)
    time_points = [current_time]
    node = Dict{Int64,Any}()
    pipe = Dict{Int64,Any}()
    compressor = Dict{Int64,Any}()

    if params(ts, :load_adjust) == true
        for i in keys(get(ref(ts), :node, []))
            node[i] = Dict(
                "pressure" => [ref(ts, :node, i, "pressure")],
                "load_reduction" => [ref(ts, :node, i, "load_reduction")],
            )
        end
    else
        for i in keys(get(ref(ts), :node, []))
            node[i] = Dict("pressure" => [ref(ts, :node, i, "pressure")])
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
        compressor[i] = Dict("flow" => [ref(ts, :compressor, i, "flow")])
    end
    return OutputState(time_points, node, pipe, compressor, current_time + output_dt)
end

function capture_step_state(ts::TransientSimulator)::StepState
    pressure_time = ref(ts, :current_time)
    flux_time = max(params(ts, :t_0), pressure_time - params(ts, :dt)/2)
    node_pressure = Dict{Int64,Float64}()
    node_load_reduction = Dict{Int64,Float64}()
    pipe_fr_mass_flux = Dict{Int64,Float64}()
    pipe_to_mass_flux = Dict{Int64,Float64}()
    compressor_flow = Dict{Int64,Float64}()

    for i in keys(get(ref(ts), :node, []))
        node_pressure[i] = ref(ts, :node, i, "pressure")
        node_load_reduction[i] = ref(ts, :node, i, "load_reduction")
    end

    for i in keys(get(ref(ts), :pipe, []))
        pipe_fr_mass_flux[i] = ref(ts, :pipe, i, "fr_mass_flux")
        pipe_to_mass_flux[i] = ref(ts, :pipe, i, "to_mass_flux")
    end

    for i in keys(get(ref(ts), :compressor, []))
        compressor_flow[i] = ref(ts, :compressor, i, "flow")
    end

    return StepState(
        pressure_time,
        flux_time,
        node_pressure,
        node_load_reduction,
        pipe_fr_mass_flux,
        pipe_to_mass_flux,
        compressor_flow,
    )
end

function _linear_interpolate(x::Float64, x0::Float64, y0::Float64, x1::Float64, y1::Float64)
    if isapprox(x0, x1; atol = TOL, rtol = TOL)
        return y1
    end
    return y0 + (y1 - y0) * ((x - x0) / (x1 - x0))
end

function _append_output_checkpoint!(
    ts::TransientSimulator,
    state::OutputState,
    output_time::Float64,
    previous_step::StepState,
    current_step::StepState,
)
    push!(state.time_points, output_time)
    for i in keys(get(ref(ts), :node, []))
        push!(
            state.node[i]["pressure"],
            _linear_interpolate(
                output_time,
                previous_step.pressure_time,
                previous_step.node_pressure[i],
                current_step.pressure_time,
                current_step.node_pressure[i],
            ),
        )
    end

    if params(ts, :load_adjust) == true
        for i in keys(get(ref(ts), :node, []))
            push!(
                state.node[i]["load_reduction"],
                _linear_interpolate(
                    output_time,
                    previous_step.flux_time,
                    previous_step.node_load_reduction[i],
                    current_step.flux_time,
                    current_step.node_load_reduction[i],
                ),
            )
        end
    end

    for i in keys(get(ref(ts), :pipe, []))
        push!(
            state.pipe[i]["fr_mass_flux"],
            _linear_interpolate(
                output_time,
                previous_step.flux_time,
                previous_step.pipe_fr_mass_flux[i],
                current_step.flux_time,
                current_step.pipe_fr_mass_flux[i],
            ),
        )
        push!(
            state.pipe[i]["to_mass_flux"],
            _linear_interpolate(
                output_time,
                previous_step.flux_time,
                previous_step.pipe_to_mass_flux[i],
                current_step.flux_time,
                current_step.pipe_to_mass_flux[i],
            ),
        )
    end
    for i in keys(get(ref(ts), :compressor, []))
        push!(
            state.compressor[i]["flow"],
            _linear_interpolate(
                output_time,
                previous_step.flux_time,
                previous_step.compressor_flow[i],
                current_step.flux_time,
                current_step.compressor_flow[i],
            ),
        )
    end
    return
end

function update_output_state!(
    ts::TransientSimulator,
    state::OutputState,
    previous_step::StepState,
    current_step::StepState;
    finalize::Bool = false,
)
    upper_time =
        finalize == true ? min(current_step.pressure_time, params(ts, :t_f)) : current_step.flux_time

    while state.next_output_time <= upper_time + TOL
        _append_output_checkpoint!(ts, state, state.next_output_time, previous_step, current_step)
        state.next_output_time += params(ts, :output_dt)
    end
    return
end

function update_output_data!(ts::TransientSimulator, state::OutputState, data::OutputData)
    data.initial_time = params(ts, :t_0)
    data.final_time = ref(ts, :current_time)
    data.time_points = copy(state.time_points)
    for (i, _) in ref(ts, :node)
        data.node[i]["pressure"] = copy(state.node[i]["pressure"])
        data.final_state["nodal_pressure"][i] = ref(ts, :node, i, "pressure")
    end

    if params(ts, :load_adjust) == true
        for (i, _) in ref(ts, :node)
            data.node[i]["load_reduction"] = copy(state.node[i]["load_reduction"])
        end
    end

    for (i, pipe) in ref(ts, :pipe)
        data.pipe[i]["fr_mass_flux"] = copy(state.pipe[i]["fr_mass_flux"])
        data.pipe[i]["to_mass_flux"] = copy(state.pipe[i]["to_mass_flux"])

        n = pipe["num_discretization_points"]
        dx = pipe["dx"]
        L = pipe["length"]
        area = pipe["area"]
        x_rho = LinRange(0, L, n)
        if params(ts, :method) in [:explicit_staggered_grid, :explicit_staggered_grid_new]
            x_mid = x_rho[1:(n-1)] .+ dx/2.0
            # this is what needs to be done to replicate
            # x_phi = [-(dx/2), x_mid..., L+(dx/2)] 
            x_phi = [0.0, x_mid..., L]
        else
            x_phi = x_rho
        end
        rho = pipe["density_profile"] # len n
        pressure = [get_pressure(ts, val) for val in rho]
        phi = pipe["mass_flux_profile"] # len n+1
        flow = phi .* area

        data.final_state["pipe_flow"][i] = Spline1D(x_phi, flow, k = 1)
        data.final_state["pipe_pressure"][i] = Spline1D(x_rho, pressure, k = 1)
    end
    for i in keys(get(ref(ts), :compressor, []))
        data.compressor[i]["flow"] = copy(state.compressor[i]["flow"])
        data.final_state["compressor_flow"][i] = ref(ts, :compressor, i, "flow")
    end
end

function update_output_data_final_state_only!(ts::TransientSimulator, data::OutputData)
    data.initial_time = params(ts, :t_0)
    data.final_time = ref(ts, :current_time)
    for (i, _) in ref(ts, :node)
        data.final_state["nodal_pressure"][i] = ref(ts, :node, i, "pressure")
    end


    for (i, pipe) in ref(ts, :pipe)
        n = pipe["num_discretization_points"]
        dx = pipe["dx"]
        L = pipe["length"]
        area = pipe["area"]
        x_rho = LinRange(0, L, n)
        if params(ts, :method) == :explicit_staggered_grid
            x_mid = x_rho[1:(n-1)] .+ dx/2.0
            # this is what needs to be done to replicate
            # x_phi = [-(dx/2), x_mid..., L+(dx/2)] 
            x_phi = [0.0, x_mid..., L]
        else
            x_phi = x_rho
        end
        rho = pipe["density_profile"] # len n
        pressure = [get_pressure(ts, val) for val in rho]
        phi = pipe["mass_flux_profile"] # len n+1
        flow = phi .* area

        data.final_state["pipe_flow"][i] = Spline1D(x_phi, flow, k = 1)
        data.final_state["pipe_pressure"][i] = Spline1D(x_rho, pressure, k = 1)
    end
    for i in keys(get(ref(ts), :compressor, []))
        data.final_state["compressor_flow"][i] = ref(ts, :compressor, i, "flow")
    end
end

function populate_solution!(ts::TransientSimulator, output::OutputData)
    times = output.time_points
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
        key = isa(i, String) ? parse(Int64, i) : i
        sol["nodes"][i]["pressure"] = pressure_convertor.(output.node[key]["pressure"])
        sol["final_state"]["initial_nodal_pressure"][i] =
            pressure_convertor(output.final_state["nodal_pressure"][key])
    end

    if params(ts, :load_adjust) == true
        for (i, _) in ref(ts, :node)
            sol["nodes"][string(i)]["load_reduction"] = Vector{Float64}()
            sol["nodes"][string(i)]["load_reduction"] =
                flow_convertor.(output.node[i]["load_reduction"])
        end
        sol["load_reduction_nodes"] = Vector{String}()
        for node_id in ref(ts, :load_reduction_nodes)
            push!(sol["load_reduction_nodes"], string(node_id))
        end
    end

    sol["final_state"]["initial_pipe_flow"] = Dict{String,Any}()
    sol["final_state"]["initial_pipe_pressure"] = Dict{String,Any}()

    for (i, pipe) in get(data, "pipes", [])
        dx = params(ts, :output_dx)
        pipe_length = pipe["length"]
        num_discretization_points = Int(floor(pipe_length/dx)) + 1
        distance = collect(range(0, pipe_length, length = num_discretization_points))
        key = isa(i, String) ? parse(Int64, i) : i
        area = ref(ts, :pipe, key, "area")
        fr_flow = output.pipe[key]["fr_mass_flux"] .* area
        to_flow = output.pipe[key]["to_mass_flux"] .* area
        sol["pipes"][i]["in_flow"] = flow_convertor.(fr_flow)
        sol["pipes"][i]["out_flow"] = flow_convertor.(to_flow)
        fr_node = string(ref(ts, :pipe, key, "fr_node"))
        to_node = string(ref(ts, :pipe, key, "to_node"))
        sol["pipes"][i]["in_pressure"] = sol["nodes"][fr_node]["pressure"]
        sol["pipes"][i]["out_pressure"] = sol["nodes"][to_node]["pressure"]
        flow_spl = output.final_state["pipe_flow"][key]
        flow = [flow_spl(x) for x in distance]
        pressure_spl = output.final_state["pipe_pressure"][key]
        pressure = [pressure_spl(x) for x in distance]
        sol["final_state"]["initial_pipe_flow"][i] = Dict{String,Any}(
            "distance" => length_convertor.(distance),
            "value" => flow_convertor.(flow),
        )
        sol["final_state"]["initial_pipe_pressure"][i] = Dict{String,Any}(
            "distance" => length_convertor.(distance),
            "value" => pressure_convertor.(pressure),
        )
    end
    sol["final_state"]["initial_compressor_flow"] = Dict{String,Any}()
    for (i, _) in get(data, "compressors", [])
        key = isa(i, String) ? parse(Int64, i) : i
        fr_node = string(ref(ts, :compressor, key, "fr_node"))
        to_node = string(ref(ts, :compressor, key, "to_node"))
        sol["compressors"][i]["suction_pressure"] = sol["nodes"][fr_node]["pressure"]
        sol["compressors"][i]["discharge_pressure"] = sol["nodes"][to_node]["pressure"]
        sol["compressors"][i]["compression_ratio"] =
            sol["nodes"][to_node]["pressure"] ./ sol["nodes"][fr_node]["pressure"]
        sol["compressors"][i]["flow"] = flow_convertor.(output.compressor[key]["flow"])
        sol["final_state"]["initial_compressor_flow"][i] =
            flow_convertor(output.final_state["compressor_flow"][key])
    end

    sol["base_time_step"] = time_convertor(params(ts, :base_dt))
    sol["time_points"] = time_convertor.(times)
    return
end

function populate_solution_final_state_only!(ts::TransientSimulator, output::OutputData)
    times = output.time_points
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
        key = isa(i, String) ? parse(Int64, i) : i
        sol["final_state"]["initial_nodal_pressure"][i] =
            pressure_convertor(output.final_state["nodal_pressure"][key])
    end



    sol["final_state"]["initial_pipe_flow"] = Dict{String,Any}()
    sol["final_state"]["initial_pipe_pressure"] = Dict{String,Any}()

    for (i, pipe) in get(data, "pipes", [])
        dx = params(ts, :output_dx)
        pipe_length = pipe["length"]
        num_discretization_points = Int(floor(pipe_length/dx)) + 1
        distance = collect(range(0, pipe_length, length = num_discretization_points))
        key = isa(i, String) ? parse(Int64, i) : i
        flow_spl = output.final_state["pipe_flow"][key]
        flow = [flow_spl(x) for x in distance]
        pressure_spl = output.final_state["pipe_pressure"][key]
        pressure = [pressure_spl(x) for x in distance]
        sol["final_state"]["initial_pipe_flow"][i] = Dict{String,Any}(
            "distance" => length_convertor.(distance),
            "value" => flow_convertor.(flow),
        )
        sol["final_state"]["initial_pipe_pressure"][i] = Dict{String,Any}(
            "distance" => length_convertor.(distance),
            "value" => pressure_convertor.(pressure),
        )
    end
    sol["final_state"]["initial_compressor_flow"] = Dict{String,Any}()
    for (i, _) in get(data, "compressors", [])
        key = isa(i, String) ? parse(Int64, i) : i
        sol["final_state"]["initial_compressor_flow"][i] =
            flow_convertor(output.final_state["compressor_flow"][key])
    end

    sol["base_time_step"] = time_convertor(params(ts, :base_dt))
    sol["time_points"] = time_convertor.(times)
    return
end
