# GasTranSim.jl Documentation

```@meta
CurrentModule = GasTranSim
```

## Calling the simulator

We shall illustrate the simplest call  to invoke the simulator and then explain the various input arguments with reference to the simple case of  a single pipe. The input files are contained in the directory `examples/data/model1pipe_fast_transients/`


The initialization function for the simulator has the definition
`initialize_simulator(data_folder::AbstractString;case_name::AbstractString="", case_types::Vector{Symbol}=Symbol[], kwargs...)`

From Julia REPL, the simplest call is:
```julia
folder = "examples/data/model1pipe_fast_transients/"
ts = initialize_simulator(folder)
# equivalent to the call
ts = initialize_simulator(folder; eos = :ideal, case_name="", case_types=[:params, :network, :bc, :ic])
```
The minimal function call implicitly looks for all input files in the path `folder`. It assumes the equation of the state for the natural gas to be that given by the ideal gas equation, and that the `case_name` is an empty string. Thus it looks for input files `params.json`, `network.json`, `bc.json`, `ic.json`.


However, suppose one wanted to use the same network, parameters and initial conditions with different boundary conditions in say the file `bc_new_case.json`, and with the CNGA equation of state. Then one would use the call:
```julia
folder = "examples/data/model1pipe_fast_transients/"
ts = initialize_simulator(folder; eos = :simple_cnga, case_name="new_case", case_types=[:bc])
```
If the initial conditions are also different, and stored in `ic_new_case.json`, then the call would be modified to
```julia
folder = "examples/data/model1pipe_fast_transients/"
ts = initialize_simulator(folder; eos = :simple_cnga, case_name="new_case", case_types=[:bc, :ic])
```

Once the simulator is initialized, it is run by invoking

```julia
run_simulator!(ts)
```


## Saving output
At the end of the simulation, there are two types of output we are interested in.

1. A record of pressures/density at every node (junction) of the network, and the inflows/outflows into each pipe of the network at a frequency (*Output dt*) specified in the input file `params.json`.
2. The state of the network at the end of the simulation, i.e., the values of the pressure/density and the mass flux along the pipes throughout the network are recorded following the spatial resolution (*Output dx*) set in the input file `params.json`. This information is stored in the same format as `ic.json` so that a new simulation can use the final state as the new initial condition and move forward in time.



Both of these tasks are accomplished as follows:

```julia
folder = "examples/data/model1pipe_fast_transients/"
ts = initialize_simulator(folder)
run_simulator!(ts)
write_output(ts; output_path = folder, 
	output_file = "output_time_history.json", final_state_file = "ic_restart.json")
```
Note that the flag for *Save final state* in `params.json` must be set to 1 in order to enable saving the output.

The file `test/sim_restart.jl` demonstrates that one run of the simulator from 0 to 24 hours is equivalent to  two successive runs, one from 0 to 12 hours, and then the next one from 12 - 24 hours.


