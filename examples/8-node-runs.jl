using GasTranSim

base_path = split(Base.active_project(), "Project.toml")[1]
folder = base_path * "data/8-node/"

ts = initialize_simulator(folder)
run_simulator!(ts; showprogress = true)
@info "First run completed"

t = ts.sol["time_points"] / 3600.0 # hrs
dt1 = params(ts, :dt) * nominal_values(ts, :time)

cratio_1 = ts.sol["compressors"]["1"]["compression_ratio"]
cratio_2 = ts.sol["compressors"]["2"]["compression_ratio"]
cratio_3 = ts.sol["compressors"]["3"]["compression_ratio"]

pressure_node_5 = ts.sol["nodes"]["5"]["pressure"]

in_flow_node_6 = ts.sol["pipes"]["1"]["in_flow"]
out_flow_node_2 = ts.sol["pipes"]["1"]["out_flow"]

out_flow_node_5 = ts.sol["pipes"]["5"]["out_flow"]
out_flow_node_3 = ts.sol["pipes"]["2"]["out_flow"] .- ts.sol["pipes"]["3"]["in_flow"]

# run Ideal Eos case 
ts = initialize_simulator(
    folder;
    eos = :ideal,
    case_name = "fine_time_step",
    case_types = [:params],
)
run_simulator!(ts)

t1 = ts.sol["time_points"] / 3600.0 # hrs
dt2 = params(ts, :dt) * nominal_values(ts, :time)

pressure_node_5_fine = ts.sol["nodes"]["5"]["pressure"]

in_flow_node_6_fine = ts.sol["pipes"]["1"]["in_flow"]
out_flow_node_2_fine = ts.sol["pipes"]["1"]["out_flow"]
@info "Second run complete"