@testset "8-node steady BC" begin
    folder = "./data/8-node-steady/"
    ts = initialize_simulator(folder)
    run_simulator!(ts)
    sol = ts.sol
    for i in keys(get(sol, "nodes", []))
        pressure = sol["nodes"][i]["pressure"]
        error = maximum( abs.(pressure .- pressure[1]) ) / pressure[1]
        @test error ≈ 0.0 atol=1e-2
    end
end 




