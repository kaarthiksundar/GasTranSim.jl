function initialize_solution(data::Dict{String,Any}, params::Dict{Symbol,Any})::Dict{String,Any}
    sol = Dict{String,Any}()
    sol["time_points"] = Vector{Float64}()
    sol["is_si_units"] = 0
    sol["is_english_units"] = 0
    sol["nodes"] = Dict{String,Any}()
    sol["pipes"] = Dict{String,Any}()
    sol["compressors"] = Dict{String,Any}()

    for (i, node) in get(data, "nodes", [])
        sol["nodes"][i] = Dict{String,Any}()
        sol["nodes"][i]["injection"] = Vector{Float64}()
        sol["nodes"][i]["pressure"] = Vector{Float64}()
    end

    for (i, pipe) in get(data, "pipes", [])
        sol["pipes"][i] = Dict{String,Any}()
        sol["pipes"][i]["in_mass_flow"] = Vector{Float64}()
        sol["pipes"][i]["out_mass_flow"] = Vector{Float64}()
        sol["pipes"][i]["in_pressure"] = Vector{Float64}()
        sol["pipes"][i]["out_pressure"] = Vector{Float64}()
    end

    for (i, compressor) in get(data, "compressors", [])
        sol["compressors"][i] = Dict{String,Any}()
        sol["compressors"][i]["mass_flow"] = Vector{Float64}()
        sol["compressors"][i]["suction_pressure"] = Vector{Float64}()
        sol["compressors"][i]["discharge_pressure"] = Vector{Float64}()
    end

    return sol
end
