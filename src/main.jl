using JSON 
using Dierckx

include("io/json.jl")
include("io/data_utils.jl")

file = "./data/model8ts_3d.json";

data = parse_json(file)

params, nominal_values = process_data!(data)