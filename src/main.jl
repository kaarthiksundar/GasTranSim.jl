using JSON
using Dierckx
using ProgressMeter

include("io/json.jl")
include("io/data_utils.jl")

include("unit_conversion/unit_convertor_utils.jl")
include("unit_conversion/to_si.jl")
include("unit_conversion/to_english.jl")
include("unit_conversion/to_pu.jl")
include("unit_conversion/unit_convertors.jl")

include("core/types.jl")
include("core/ref.jl")
include("core/bc.jl")
include("core/sol.jl")
include("core/initialize_ts.jl")
include("core/run_task.jl")
include("core/time_integration.jl")
include("core/run_ts.jl")

include("io/output.jl")

# file = "./data/model8ts_3d.json";
file = "./data/model30_ts.json"

ts = initialize_simulator(file)

run_simulator!(ts)
