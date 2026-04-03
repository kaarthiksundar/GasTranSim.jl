@testset "Output checkpointing logic" begin
    invalid_ts =
        initialize_simulator("./data/8-node/"; case_name = "first_half", case_types = [:params])
    invalid_ts.params[:output_dt] = 0.0
    @test_throws DomainError GasTranSim.run_simulator!(
        invalid_ts;
        showprogress = false,
        turnoffprogressbar = true,
    )

    ts =
        initialize_simulator("./data/8-node/"; case_name = "first_half", case_types = [:params])
    ts.params[:load_adjust] = true
    ts.params[:t_f] = 2.0
    ts.params[:output_dt] = 1.0

    state = GasTranSim.initialize_output_state(ts)
    initial_pressure = state.node[1]["pressure"][1]
    initial_load_reduction = state.node[1]["load_reduction"][1]
    initial_fr_mass_flux = state.pipe[1]["fr_mass_flux"][1]
    initial_to_mass_flux = state.pipe[1]["to_mass_flux"][1]
    initial_compressor_flow = state.compressor[1]["flow"][1]

    node_ids = collect(keys(GasTranSim.ref(ts, :node)))
    pipe_ids = collect(keys(GasTranSim.ref(ts, :pipe)))
    compressor_ids = collect(keys(get(GasTranSim.ref(ts), :compressor, Dict{Int64,Any}())))

    previous_step = GasTranSim.StepState(
        0.0,
        0.0,
        Dict(i => 10.0 for i in node_ids),
        Dict(i => 1.0 for i in node_ids),
        Dict(i => 2.0 for i in pipe_ids),
        Dict(i => 4.0 for i in pipe_ids),
        Dict(i => 6.0 for i in compressor_ids),
    )

    current_step = GasTranSim.StepState(
        2.0,
        1.5,
        Dict(i => 14.0 for i in node_ids),
        Dict(i => 3.0 for i in node_ids),
        Dict(i => 6.0 for i in pipe_ids),
        Dict(i => 8.0 for i in pipe_ids),
        Dict(i => 10.0 for i in compressor_ids),
    )

    GasTranSim.update_output_state!(ts, state, previous_step, current_step)

    @test state.time_points == [0.0, 1.0]
    @test state.node[1]["pressure"][1] == initial_pressure
    @test state.node[1]["pressure"][2] ≈ 12.0
    @test state.node[1]["load_reduction"][1] == initial_load_reduction
    @test state.node[1]["load_reduction"][2] ≈ 7.0 / 3.0
    @test state.pipe[1]["fr_mass_flux"][1] == initial_fr_mass_flux
    @test state.pipe[1]["fr_mass_flux"][2] ≈ 14.0 / 3.0
    @test state.pipe[1]["to_mass_flux"][1] == initial_to_mass_flux
    @test state.pipe[1]["to_mass_flux"][2] ≈ 20.0 / 3.0
    @test state.compressor[1]["flow"][1] == initial_compressor_flow
    @test state.compressor[1]["flow"][2] ≈ 26.0 / 3.0
    @test state.next_output_time == 2.0
end
