using GasTranSim
# import PyPlot
# const plt = PyPlot
# using LaTeXStrings
# using TerminalExtensions




folder = "./data/8-node-ramp-steady/"
##=== Run from t=0 to t=t1 ===##

ts = initialize_simulator(folder; case_name="ramp", case_types=[:bc, :ic]);
run_simulator!(ts, turnoffprogressbar=true)
# write_output(ts; output_path = folder)


println("run complete")




