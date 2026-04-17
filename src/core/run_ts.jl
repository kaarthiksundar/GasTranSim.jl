function run_simulator!(
    ts::TransientSimulator;
    run_type::Symbol = :serial,
    save_snapshots::Bool = false,
    snapshot_percent::Float64 = 10.0,
    snapshot_path::AbstractString = "./",
    snapshot_filename::AbstractString = "solution-snapshot",
    steady_state::Bool = false,
    load_adjust::Bool = false,
    showprogress::Bool = true,
    turnoffprogressbar::Bool = false,
    progress_dt = 1.0,
)
    minimum_pressure_limit = params(ts, :minimum_pressure_limit)
    output_dt = params(ts, :output_dt)
    base_dt = params(ts, :base_dt)
    ts.params[:load_adjust] = load_adjust

    if output_dt <= 0.0
        throw(DomainError(output_dt, "output_dt must be > 0"))
    end

    if base_dt <= 0.0
        throw(DomainError(base_dt, "base_dt must be > 0"))
    end

    if params(ts, :t_f) < params(ts, :t_0)
        throw(DomainError(params(ts, :t_f), "t_f must be >= t_0"))
    end

    if save_snapshots == true && snapshot_percent <= 0.0
        throw(
            DomainError(snapshot_percent, "snapshot_percent must be > 0 when saving snapshots"),
        )
    end

    if params(ts, :load_adjust) == true && !(minimum_pressure_limit > 0)
        throw(
            DomainError(
                minimum_pressure_limit,
                "load adjustment requires minimum_pressure_limit > 0",
            ),
        )
    end
    #
    (params(ts, :load_adjust) == true) && (ts.ref[:load_reduction_nodes] = Vector{Int64}())

    output_state = initialize_output_state(ts)
    t_f = params(ts, :t_f)
    t_0 = params(ts, :t_0)
    total_time = t_f - t_0
    #
    output_data = OutputData(ts)
    previous_step_state = capture_step_state(ts)
    #
    progress_total = 1000
    prog = Progress(
        progress_total;
        dt = progress_dt,
        barglyphs = get_barglyphs(),
        barlen = 10,
        enabled = showprogress,
        desc = "Sim. progress: ",
    )
    if showprogress == false
        prog = ProgressUnknown(desc = "Sim. status ", spinner = true)
    end
    # This block is used only for computing steady-state solution
    if steady_state == true
        nodal_pressure_previous = form_nodal_pressure_vector(ts)
        @info "Change in nodal pressure will be computed after 10%, 20%,...100% of total time"
    end
    # Saving snapshot of initial condition
    snapshot_count = 0
    snapshot_interval = (t_f - t_0) * snapshot_percent / 100.0
    next_snapshot_time = t_0 + snapshot_interval
    steady_state_check_interval = total_time / 10.0
    next_steady_state_check_time = t_0 + steady_state_check_interval
    if save_snapshots == true
        snapshot_count =
            save_snapshot(ts, output_data, snapshot_path, snapshot_filename, snapshot_count)
    end

    lin_system = form_matrix_for_compressor_flow_solve(ts)
    # Time marching loop
    step = 0
    while ref(ts, :current_time) < t_f - TOL
        step += 1
        dt_step = min(base_dt, t_f - ref(ts, :current_time))
        ts.params[:dt] = dt_step
        #
        advance_current_time!(ts, dt_step)
        #  if current_time is where some disruption occurs, modify ts.ref now
        advance_pipe_density_internal!(ts, run_type) # (n+1) level
        advance_node_pressure_mass_flux!(ts, run_type) # pressure (n+1), flux (n+1/2)
        advance_pipe_mass_flux_internal!(ts, run_type) # (n + 1 + 1/2) level
        # _compute_compressor_flows!(ts)
        _solve_compressor_flows!(ts, lin_system)
        calculate_slack_injections!(ts)
        #  if current_time is one where output needs to be saved, check and do now
        current_step_state = capture_step_state(ts)
        update_output_state!(
            ts,
            output_state,
            previous_step_state,
            current_step_state;
            finalize = (ref(ts, :current_time) >= t_f - TOL),
        )
        #
        #  This block is used only for saving snapshot of solution
        should_save_snapshot =
            save_snapshots == true &&
            (
                _should_save_snapshot_at_time(ref(ts, :current_time), t_f, next_snapshot_time) ||
                ref(ts, :current_time) >= t_f - TOL
            )
        if should_save_snapshot
            snapshot_count = save_snapshot(
                ts,
                output_data,
                snapshot_path,
                snapshot_filename,
                snapshot_count,
            )
            next_snapshot_time = _next_snapshot_time(
                ref(ts, :current_time),
                next_snapshot_time,
                snapshot_interval,
            )
        end
        previous_step_state = current_step_state
        #
        if showprogress == false
            (turnoffprogressbar == false) && (next!(prog, spinner = "🌑🌒🌓🌔🌕🌖🌗🌘"))
        else
            progress_fraction =
                (total_time <= TOL) ? 1.0 : (ref(ts, :current_time) - t_0) / total_time
            progress_fraction = min(max(progress_fraction, 0.0), 1.0)
            progress_value = round(Int, progress_fraction * progress_total)
            (turnoffprogressbar == false) && (update!(prog, progress_value))
        end
        # This block is used only for computing steady-state solution
        if steady_state == true && steady_state_check_interval > TOL
            while next_steady_state_check_time <= ref(ts, :current_time) + TOL
                nodal_pressure_current = form_nodal_pressure_vector(ts)
                error = maximum(abs.(nodal_pressure_current - nodal_pressure_previous))
                nodal_pressure_previous = nodal_pressure_current
                @info "Max change in nodal pressure: $(round(error; digits=8))"
                if error < 1e-5
                    @info "Steady state attained"
                    break
                end
                next_steady_state_check_time += steady_state_check_interval
            end
            #
            if ref(ts, :current_time) >= t_f - TOL
                @info "Steady-state not yet attained! Consider increasing final time."
            end
        end

    end
    ts.params[:dt] = base_dt
    (turnoffprogressbar == false) && (finish!(prog))
    update_output_data!(ts, output_state, output_data)
    populate_solution!(ts, output_data)

    if save_snapshots == true
        @info "Number of solution snapshots is $(snapshot_count)"
    end
end

function _should_save_snapshot_at_time(
    current_time::Float64,
    final_time::Float64,
    next_snapshot_time::Float64,
)::Bool
    if current_time >= final_time
        return false
    end
    return current_time >= next_snapshot_time
end

function _next_snapshot_time(
    current_time::Float64,
    next_snapshot_time::Float64,
    snapshot_interval::Float64,
)::Float64
    while next_snapshot_time <= current_time
        next_snapshot_time += snapshot_interval
    end
    return next_snapshot_time
end

function advance_current_time!(ts::TransientSimulator, tau::Real)
    ts.ref[:current_time] += tau
    return
end

function advance_pipe_density_internal!(ts::TransientSimulator, run_type::Symbol)
    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_advance_pipe_density_internal!, ts, key_array, run_type)
    return
end

function advance_node_pressure_mass_flux!(ts::TransientSimulator, run_type::Symbol)
    t = ref(ts, :current_time)

    # initialize a Jacobian matrix of correct size
    num_nodes = length(ref(ts, :node))
    J = spzeros(num_nodes, num_nodes)
    rhs = zeros(num_nodes)
    assemble_compressor_equations!(ts, J, rhs)
    assemble_nodal_equations!(ts, J, rhs)
    # solve nodal density at t+dt using assembled Jacobian and rhs
    x = J \ rhs
    # update nodal pressures in ts.ref using solution x
    for node_id = 1 : length(x)
        p_val = get_pressure(ts, x[node_id])
        ref(ts, :node, node_id)["pressure_previous"] = ref(ts, :node, node_id)["pressure"]
        ref(ts, :node, node_id)["pressure"] = p_val
        ref(ts, :node, node_id)["is_updated"] = true
    end

    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_compute_pipe_end_fluxes_densities!, ts, key_array, run_type)
    return
end

function advance_pipe_mass_flux_internal!(ts::TransientSimulator, run_type::Symbol)
    key_array = collect(keys(ref(ts, :pipe)))
    _execute_task!(_advance_pipe_mass_flux_internal!, ts, key_array, run_type)
    return
end

function form_nodal_pressure_vector(ts::TransientSimulator)::Vector{Float64}
    key_array = sort(collect(keys(ref(ts, :node))))
    pressure_vector = Vector{Float64}(undef, length(key_array))
    for i in eachindex(key_array)
        pressure_vector[i] = ref(ts, :node, key_array[i])["pressure"]
    end
    return pressure_vector
end

function save_snapshot(
    ts::TransientSimulator,
    output_data::OutputData,
    snapshot_path::AbstractString,
    snapshot_filename::AbstractString,
    snapshot_count::Int64,
)::Int64
    update_output_data_final_state_only!(ts, output_data)
    populate_solution_final_state_only!(ts, output_data)
    write_final_state(
        ts;
        output_path = snapshot_path,
        final_state_file = "$(snapshot_filename)-$(snapshot_count).json",
    )
    return snapshot_count + 1
end


