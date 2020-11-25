struct TransientSimulator 
    data::Dict{String,Any}
    ref::Dict{Symbol,Any}
    sol::Dict{String,Any}
    nominal_values::Dict{Symbol,Any}
    params::Dict{Symbol,Any}
    boundary_conditions::Dict{Symbol,Any}
end 

ref(ts::TransientSimulator) = ts.ref 
ref(ts::TransientSimulator, key::Symbol) = ts.ref[key]
ref(ts::TransientSimulator, key::Symbol, id::Int64) = ts.ref[key][id]
ref(ts::TransientSimulator, key::Symbol, id::Int64, field) = ts.ref[key][id][field]

params(ts::TransientSimulator) = ts.params 
params(ts::TransientSimulator, key::Symbol) = ts.params[key]

nominal_values(ts::TransientSimulator) = ts.params 
nominal_values(ts::TransientSimulator, key::Symbol) = ts.params[key]

function control(ts::TransientSimulator, 
    key::Symbol, id::Int64, t::Real)::Tuple{CONTROL_TYPE,Float64}
    (key == :node) && (return get_nodal_control(ts, id, t))
    (key == :compressor) && (return get_compressor_control(ts, id, t))
    @error "control available only for nodes and compressors"
    return CONTROL_TYPE::unknown, 0.0
end 

function get_nodal_control(ts::TransientSimulator, 
    id::Int64, t::Real)::Tuple{CONTROL_TYPE,Float64}
    if !haskey(ts.boundary_conditions[:node], id)
        return CONTROL_TYPE::flow, 0.0
    end 
    val = ts.boundary_conditions[:node][id]["spl"](t)
    control_type = ts.boundary_conditions[:node][id]["control_type"]
    return control_type, val
end 

function get_compressor_control(ts::TransientSimulator, 
    id::Int64, t::Real)::Tuple{CONTROL_TYPE,Float64}
    control_type, val = ts.boundary_conditions[:compressor][id]["spl"](t)
    return CONTROL_TYPE(control_type), val
end 

struct CompressorControl 
    time::Vector{Float64}
    control_type::Vector{Any}
    value::Vector{Float64}
end 

function evaluate(sp::CompressorControl, t::Real)::Tuple{Any,Float64}
    @assert (t >= sp.time[1]) && (t <= sp.time[end])
    time = sp.time 
    control_type = sp.control_type 
    value = sp.value 
    (t == time[1]) && (return control_type[1], value[1])
    (t == time[end]) && (return control_type[end], value[end])
    index = findfirst(x -> t >= x, time)
    (t == time[index]) && (return control_type[index], value[index])
    first_control = control_type[index]
    second_control = control_type[index + 1] 
    (first_control != second_control) && (return first_control, value[index])
    first_value = value[index]
    second_value = value[index + 1]
    delta_value = second_value - first_value 
    first_time = time[index]
    second_time = time[index + 1]
    delta_t = second_time - first_time  
    if (first_control == second_control)
        value = delta_value / delta_t * (t - first_time) + first_value 
        return first_control, value 
    end 
    return NaN, NaN
end 

# call synomym for CompressorControl evaluation
(spl::CompressorControl)(t::Real) = evaluate(spl, t)


@enum CONTROL_TYPE begin 
    c_ratio = 0 
    discharge_pressure = 1 
    flow = 2
    pressure = 10
    unknown = 100 
end 