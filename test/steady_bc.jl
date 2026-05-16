methods =[:implicit_hyperbolic, :explicit_hyperbolic, :explicit_staggered_grid, :explicit_staggered_grid_new,:implicit_parabolic]



for method in methods
    pipe_segs = nothing
    if method == :implicit_hyperbolic
        pipe_segs = 150
    end
    if method == :implicit_parabolic
        pipe_segs = 50
    end

    @info("Testing method $method...\n") 
    @testset "8-node steady BC" begin
        folder = "./data/8-node-steady/"
        ts = initialize_simulator(folder; method=method, pipe_segments=pipe_segs)
        if method == :implicit_hyperbolic
            ts.params[:base_dt] = 100 * ts.params[:base_dt]
        end
        if method == :implicit_parabolic
            ts.params[:base_dt] = 100 * ts.params[:base_dt]
        end
        run_simulator!(ts; method=method, turnoffprogressbar = true)
        sol = ts.sol
        for i in keys(get(sol, "nodes", []))
            pressure = sol["nodes"][i]["pressure"]
            error = maximum(abs.(pressure .- pressure[1])) / pressure[1]
            @test error ≈ 0.0 atol=1.5e-2
        end
        ts = initialize_simulator(folder; method=method, pipe_segments=pipe_segs, eos = :simple_cnga)
        if method == :implicit_hyperbolic
            ts.params[:base_dt] = 100 * ts.params[:base_dt]
        end
        if method == :implicit_parabolic
            ts.params[:base_dt] = 100 * ts.params[:base_dt]
        end
        run_simulator!(ts; method=method)
        sol = ts.sol
        for i in keys(get(sol, "nodes", []))
            pressure = sol["nodes"][i]["pressure"]
            error = maximum(abs.(pressure .- pressure[1])) / pressure[1]
            @test error <= 0.06
        end

        ts = initialize_simulator(folder; method=method, pipe_segments=pipe_segs, eos = :full_cnga)
        if method == :implicit_hyperbolic
            ts.params[:base_dt] = 100 * ts.params[:base_dt]
        end
        if method == :implicit_parabolic
            ts.params[:base_dt] = 100 * ts.params[:base_dt]
        end
        run_simulator!(ts; method=method)
        sol = ts.sol
        for i in keys(get(sol, "nodes", []))
            pressure = sol["nodes"][i]["pressure"]
            error = maximum(abs.(pressure .- pressure[1])) / pressure[1]
            @test error <= 0.06
        end
    end
end
