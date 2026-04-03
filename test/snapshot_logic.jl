@testset "Snapshot scheduling logic" begin
    @test_throws DomainError GasTranSim.run_simulator!(
        initialize_simulator("./data/8-node/"; case_name = "first_half", case_types = [:params]),
        save_snapshots = true,
        snapshot_percent = 0.0,
        showprogress = false,
        turnoffprogressbar = true,
    )

    let final_time = 10.0
        snapshot_interval = final_time * 33.0 / 100.0
        next_snapshot_time = snapshot_interval
        save_times = Float64[]

        for current_time in 1.0:1.0:final_time
            should_save =
                GasTranSim._should_save_snapshot_at_time(
                    current_time,
                    final_time,
                    next_snapshot_time,
                ) || current_time == final_time
            if should_save
                push!(save_times, current_time)
                next_snapshot_time = GasTranSim._next_snapshot_time(
                    current_time,
                    next_snapshot_time,
                    snapshot_interval,
                )
            end
        end

        @test save_times == [4.0, 7.0, 10.0]
    end

    let final_time = 3.0
        snapshot_interval = final_time * 1.0 / 100.0
        next_snapshot_time = snapshot_interval
        save_times = Float64[]

        for current_time in 1.0:1.0:final_time
            should_save =
                GasTranSim._should_save_snapshot_at_time(
                    current_time,
                    final_time,
                    next_snapshot_time,
                ) || current_time == final_time
            if should_save
                push!(save_times, current_time)
                next_snapshot_time = GasTranSim._next_snapshot_time(
                    current_time,
                    next_snapshot_time,
                    snapshot_interval,
                )
            end
        end

        @test save_times == [1.0, 2.0, 3.0]
    end

    let final_time = 60.0
        snapshot_interval = final_time * 1.667 / 100.0
        next_snapshot_time = snapshot_interval
        save_times = Float64[]

        for current_time in 1.0:1.0:final_time
            should_save =
                GasTranSim._should_save_snapshot_at_time(
                    current_time,
                    final_time,
                    next_snapshot_time,
                ) || current_time == final_time
            if should_save
                push!(save_times, current_time)
                next_snapshot_time = GasTranSim._next_snapshot_time(
                    current_time,
                    next_snapshot_time,
                    snapshot_interval,
                )
            end
        end

        @test save_times == collect(2.0:1.0:60.0)
    end
end
