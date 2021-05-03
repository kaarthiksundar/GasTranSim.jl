function initialize_solution(data::Dict{String,Any}, params::Dict{Symbol,Any})::Dict{String,Any}
    sol = Dict{String,Any}()
    sol["time_points"] = Vector{Float64}()
    sol["is_si_units"] = 0
    sol["is_english_units"] = 0
    sol["is_per_unit"] = 1
    sol["nodes"] = Dict{String,Any}()
    sol["pipes"] = Dict{String,Any}()
    sol["compressors"] = Dict{String,Any}()

    for (i, node) in get(data, "nodes", [])
        sol["nodes"][i] = Dict{String,Any}()
        sol["nodes"][i]["flow"] = Vector{Float64}()
        sol["nodes"][i]["pressure"] = Vector{Float64}()
    end

    for (i, pipe) in get(data, "pipes", [])
        sol["pipes"][i] = Dict{String,Any}()
        sol["pipes"][i]["in_flow"] = Vector{Float64}()
        sol["pipes"][i]["out_flow"] = Vector{Float64}()
        sol["pipes"][i]["in_pressure"] = Vector{Float64}()
        sol["pipes"][i]["out_pressure"] = Vector{Float64}()
    end

    for (i, compressor) in get(data, "compressors", [])
        sol["compressors"][i] = Dict{String,Any}()
        sol["compressors"][i]["flow"] = Vector{Float64}()
        sol["compressors"][i]["suction_pressure"] = Vector{Float64}()
        sol["compressors"][i]["discharge_pressure"] = Vector{Float64}()
    end

    return sol
end

function update_solution!(ts::TransientSimulator) 
    push!(solution["time_points"], ref(ts, :current_time))

    for (i, node) in get(data, "nodes", [])
        key = isa(i, Number) ? i : parse(Int64, i) 
        push!(sol["nodes"][i]["pressure"], ref(ts, :node, key, "pressure"))
        push!(sol["nodes"][i]["flow"], ref(ts, :node, key, "injection"))
    end 

end 