"Grab the data from a json field"
function parse_json(file_string::AbstractString)
    (~isfile(file_string)) && (return Dict{String,Any}())
    data = open(file_string, "r") do io
        parse_json(io)
    end
    return data
end


""
function parse_json(io::IO)
    data = JSON.parse(io, dicttype = Dict)
    return data
end
