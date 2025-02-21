@testset "Simulation restart" begin
    folder = "./data/8-node/"
    # full run for 24 hours
    ts = initialize_simulator(folder; 
        case_name="full", case_types=[:params])
    run_simulator!(ts, run_type=:serial, showprogress=true)
    pressure_node_5 = ts.sol["nodes"]["5"]["pressure"]
    inflow_node_6 = ts.sol["pipes"]["1"]["in_flow"]
    outflow_node_2 = ts.sol["pipes"]["1"]["out_flow"]

    # run for 0 to 12 hours 
    ts_a = initialize_simulator(folder; 
        case_name="first_half", case_types=[:params])
    run_simulator!(ts_a)
    final_state = Dict(
        key => value for (key, value) in ts_a.sol["final_state"]
    )
    pressure_node_5_a = ts_a.sol["nodes"]["5"]["pressure"]
    inflow_node_6_a = ts_a.sol["pipes"]["1"]["in_flow"]
    outflow_node_2_a = ts_a.sol["pipes"]["1"]["out_flow"]


    # run for 12 to 24 hours 
    ts_b = initialize_simulator(folder; 
        case_name="second_half", case_types=[:params, :ic])
    run_simulator!(ts_b, showprogress=false)
    final_state = Dict(
        key => value for (key, value) in ts_b.sol["final_state"]
    )
    pressure_node_5_b = ts_b.sol["nodes"]["5"]["pressure"]
    inflow_node_6_b = ts_b.sol["pipes"]["1"]["in_flow"]
    outflow_node_2_b = ts_b.sol["pipes"]["1"]["out_flow"]



    # concatenated full solution 
    pressure_node_5_cat = 
        vcat(pressure_node_5_a, pressure_node_5_b[2:end])
    inflow_node_6_cat = 
        vcat(inflow_node_6_a, inflow_node_6_b[2:end])
    outflow_node_2_cat = 
        vcat(outflow_node_2_a, outflow_node_2_b[2:end])

    # average relative errors 
    e_pressure = sum(abs.(pressure_node_5_cat - pressure_node_5) ./ 
        abs.(pressure_node_5)) / length(pressure_node_5)
    e_inflow = sum(abs.(inflow_node_6_cat - inflow_node_6) ./ 
        abs.(inflow_node_6)) / length(inflow_node_6)
    e_outflow = sum(abs.(outflow_node_2_cat - outflow_node_2) ./ 
        abs.(outflow_node_2)) / length(outflow_node_2)

    @test e_pressure ≈ 0.0 atol=1e-2
    @test e_inflow ≈ 0.0 atol=1e-2
    @test e_outflow ≈ 0.0 atol=1e-2
end 

