
using GasTranSim

# this file reproduces the runs in the paper by Gyrya, Vitaliy and Zlotnik, Anatoly

# setting up the paths for the data and outputs 
base_path = split(Base.active_project(), "Project.toml")[1]
folder = base_path * "data/1-pipe-fast-transients/"
output_json = base_path * "output/solution/"

""" 
    Note: Temp for ideal gas case set to 239.11 K, 
    nonideal 288.7K to match pressure and density
"""

# Ideal EoS run 

ts = initialize_simulator(folder; eos = :ideal)
run_simulator!(
    ts;
    save_snapshots = true,
    snapshot_path = output_json,
    snapshot_filename = "1-pipe-fast",
    snapshot_percent = 1.667,
)

nominal_density = nominal_values(ts, :density)
nominal_pressure = nominal_values(ts, :pressure)
solution = ts.sol
pipe_solution = solution["pipes"]
node_solution = solution["nodes"]

inlet_pressure = node_solution["1"]["pressure"]
inlet_density = nominal_density * get_density(ts, inlet_pressure / nominal_pressure)
outlet_pressure = node_solution["2"]["pressure"]
outlet_density = nominal_density * get_density(ts, outlet_pressure / nominal_pressure)

area = ref(ts, :pipe, 1, "area")
inlet_mass_flux = pipe_solution["1"]["in_flow"] / area
inlet_velocity = inlet_mass_flux ./ inlet_density
outlet_mass_flux = pipe_solution["1"]["out_flow"] / area
outlet_velocity = outlet_mass_flux ./ outlet_density

t = ts.sol["time_points"] / 60.0 # mins
@info "ideal run completed"

# CNGA EoS run 

ts_cnga = initialize_simulator(
    folder;
    eos = :simple_cnga,
    case_name = "cnga",
    case_types = [:params],
)
run_simulator!(ts_cnga, save_snapshots = false)

t_cnga = ts_cnga.sol["time_points"] / 60.0 # mins

nominal_density = nominal_values(ts_cnga, :density)
nominal_pressure = nominal_values(ts_cnga, :pressure)
solution = ts_cnga.sol
pipe_solution = solution["pipes"]
node_solution = solution["nodes"]

inlet_pressure_cnga = node_solution["1"]["pressure"]
inlet_density_cnga =
    nominal_density * get_density(ts_cnga, inlet_pressure_cnga / nominal_pressure)
outlet_pressure_cnga = node_solution["2"]["pressure"]
outlet_density_cnga =
    nominal_density * get_density(ts_cnga, outlet_pressure_cnga / nominal_pressure)

area = ref(ts, :pipe, 1, "area")
inlet_mass_flux_cnga = pipe_solution["1"]["in_flow"] / area
inlet_velocity_cnga = inlet_mass_flux_cnga ./ inlet_density_cnga
outlet_mass_flux_cnga = pipe_solution["1"]["out_flow"] / area
outlet_velocity_cnga = outlet_mass_flux_cnga ./ outlet_density_cnga

@info "simple CNGA run completed"
