using GasTranSim

using Test

GasTranSim.set_logging_level!(:Debug)

include("units.jl")

include("sim_restart.jl")
    
include("steady_bc.jl")

include("logging_tests.jl")
