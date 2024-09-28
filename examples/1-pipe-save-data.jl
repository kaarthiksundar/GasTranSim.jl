
using GasTranSim


base_path = split(Base.active_project(), "Project.toml")[1]
folder = base_path * "data/1-pipe-fast-transients/"
output_json = base_path * "output/solution/"

""" 
    Note: Temp for ideal gas case set to 239.11 K, 
    nonideal 288.7K to match pressure and density
"""

## === Save final state at every minute === ##
@info "runs started"
for t_f in 60:60:3600
    ts = initialize_simulator(folder; eos = :ideal)
    params(ts)[:t_f] = float(t_f) / nominal_values(ts, :time)
    run_simulator!(ts)
    t_f_min = Int(t_f/60)
    write_final_state(ts; output_path = output_json, 
        final_state_file = "1-pipe-fast-state-$(t_f_min)-min.json")
end
@info "runs completed"
