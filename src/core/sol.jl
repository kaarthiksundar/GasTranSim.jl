function initialize_solution(data::Dict{String,Any})::Dict{String,Any}
    sol = Dict{String,Any}()
    sol["time_step"] = Float64
    sol["time_points"] = Vector{Float64}()
    sol["nodes"] = Dict{String,Any}()
    sol["pipes"] = Dict{String,Any}()
    sol["compressors"] = Dict{String,Any}()
    sol["final_state"] = Dict{String,Any}()

    for (i, _) in get(data, "nodes", [])
        sol["nodes"][i] = Dict{String,Any}()
        sol["nodes"][i]["pressure"] = Vector{Float64}()
    end

    for (i, _) in get(data, "pipes", [])
        sol["pipes"][i] = Dict{String,Any}()
        sol["pipes"][i]["in_flow"] = Vector{Float64}()
        sol["pipes"][i]["out_flow"] = Vector{Float64}()
        sol["pipes"][i]["in_pressure"] = Vector{Float64}()
        sol["pipes"][i]["out_pressure"] = Vector{Float64}()
    end

    for (i, _) in get(data, "compressors", [])
        sol["compressors"][i] = Dict{String,Any}()
        sol["compressors"][i]["flow"] = Vector{Float64}()
        sol["compressors"][i]["suction_pressure"] = Vector{Float64}()
        sol["compressors"][i]["discharge_pressure"] = Vector{Float64}()
        sol["compressors"][i]["compression_ratio"] = Vector{Float64}()

    end

    return sol
end
