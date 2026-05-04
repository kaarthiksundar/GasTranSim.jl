methods =[:explicit_staggered_grid, :explicit_staggered_grid_new,:implicit_parabolic]
for method in methods
    @info("Testing method $method...\n") 
    @testset "8-node steady BC" begin
        folder = "./data/8-node-steady/"
        ts = initialize_simulator(folder; method=method)
        if method == :implicit_parabolic
            ts.params[:base_dt] = 100 * ts.params[:base_dt]
        end
        run_simulator!(ts; method=method, turnoffprogressbar = true)
        sol = ts.sol
        for i in keys(get(sol, "nodes", []))
            pressure = sol["nodes"][i]["pressure"]
            error = maximum(abs.(pressure .- pressure[1])) / pressure[1]
            @test error ≈ 0.0 atol=1e-2
        end
        ts = initialize_simulator(folder; method=method, eos = :simple_cnga)
        if method == :implicit_parabolic
            ts.params[:base_dt] = 100 * ts.params[:base_dt]
        end
        run_simulator!(ts; method=method)
        sol = ts.sol
        for i in keys(get(sol, "nodes", []))
            pressure = sol["nodes"][i]["pressure"]
            error = maximum(abs.(pressure .- pressure[1])) / pressure[1]
            @test error <= 0.05
        end

        ts = initialize_simulator(folder; method=method, eos = :full_cnga)
        if method == :implicit_parabolic
            ts.params[:base_dt] = 100 * ts.params[:base_dt]
        end
        run_simulator!(ts; method=method)
        sol = ts.sol
        for i in keys(get(sol, "nodes", []))
            pressure = sol["nodes"][i]["pressure"]
            error = maximum(abs.(pressure .- pressure[1])) / pressure[1]
            @test error <= 0.05
        end
    end
end
