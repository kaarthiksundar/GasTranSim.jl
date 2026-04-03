using GasTranSim

using Test

GasTranSim.set_logging_level!(:Info)

include("units.jl")

include("sim_restart.jl")
include("snapshot_logic.jl")
include("output_checkpointing.jl")

GasTranSim.set_logging_level!(:Debug)

include("steady_bc.jl")

include("logging_tests.jl")
