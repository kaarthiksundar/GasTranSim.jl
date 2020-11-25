using JSON 
using Dierckx

include("io/json.jl")
include("io/data_utils.jl")

include("unit_conversion/unit_convertor_utils.jl")
include("unit_conversion/to_si.jl")
include("unit_conversion/to_english.jl")
include("unit_conversion/to_pu.jl")
include("unit_conversion/unit_convertors.jl")

include("core/types.jl")
include("core/ref.jl")

file = "./data/model8ts_3d.json";
data = parse_json(file)
params, nominal_values = process_data!(data)

make_per_unit!(data, params, nominal_values)

ref = build_ref(data)