using GasModels 

data_file = "../data/model6.m"

data = parse_file(data_file)

ref = build_ref(data)