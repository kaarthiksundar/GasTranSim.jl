# GasTransSim.jl 

## Staged
- Adds `turnoffprogressbar` keyword argument to `run_simulator!` 
- Modifies tests to improve coverage 
- Adds new exception `MissingDataException` for missing data errors
- Adds a steady state option to run a transient simulator to get steady state solution

## v0.2.4 
- Allow for injection shedding if pressure goes above maximum pressure in params 
- Adds a `load_adjust` keyword argument to `run_ts` to activate this feature


## v0.2.3 
- Allow for load shedding if pressure falls below minimum pressure in params 
  
## v0.2.2 
- Bug fix in level computation of nodes 
- New exceptions: `ICException`, `NetworkException`, `ControlException`
- Adds more options for progress tracking
- Sets up parallel computation using `Threads.@threads`
- Adds more topology checks for two compressors in series
- Calculates initial compressor flow (if not provided) for such compressors

## v0.2.1
- Updates compat for LoggingExtras

## v0.2.0 
- Adds compressor flow computation 

## v0.1.1
- Adds GTS logger with Debug visibility on. 
- Adds progress bar flag

## v0.1.0
- Initial Release of the GasTransSim code
