using NGTransientSimulator


# folder = "./data/model8/"
# folder = "./data/model30/"
# folder = "./data/model8_paper_VG_AZ/"
folder = "./data/model1pipe_slow_transients/"
# folder = "./data/model1pipe_fast_transients/"

ts = initialize_simulator(folder, eos=:ideal)
run_simulator!(ts)
write_output(ts)

# output foldername, output filename, final state output filename
