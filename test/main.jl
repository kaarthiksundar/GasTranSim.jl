using GasTranSim


# folder = "test/data/model8/"
# folder = "test/data/model30/"
# folder = "test/data/model8_paper_VG_AZ/"
folder = "test/data/model1pipe_slow_transients/"
# folder = "test/data/model1pipe_fast_transients/"

ts = initialize_simulator(folder, eos=:ideal)
run_simulator!(ts)
write_output(ts; output_path = folder, 
	output_file = "output.json", final_state_file = "ic_final.json")

