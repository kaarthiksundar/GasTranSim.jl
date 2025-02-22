using GasTranSim

folder = "./data/elpaso/"

ts = initialize_simulator(folder; case_name="ramp", case_types=[:params,:bc, :ic]);
run_simulator!(ts, turnoffprogressbar=true, steady_state=true)
# write_output(ts; output_path = folder)

println("run complete")




