struct TransientSimulator
    data::Dict{String,Any}
    ref::Dict{Symbol,Any}
    sol::Dict{String,Any}
    nominal_values::Dict{Symbol,Any}
    params::Dict{Symbol,Any}
    initial_conditions::Dict{Symbol,Any}
    boundary_conditions::Dict{Symbol,Any}
    pu_eos_coeffs::Function
    pu_pressure_to_pu_density::Function
    pu_density_to_pu_pressure::Function
end

struct OutputState
    time_pressure::Vector{Float64}
    time_flux::Vector{Float64}
    node::Dict{Int64,Any}
    pipe::Dict{Int64,Any}
    compressor::Dict{Int64,Any}
end

mutable struct OutputData
    initial_time::Float64 
    final_time::Float64 
    node::Dict{Int64,Any}
    pipe::Dict{Int64,Any}
    compressor::Dict{Int64,Any}
    final_state::Dict{Any,Any}
end 

mutable struct CFLException <: Exception
    var::AbstractString
end

mutable struct NetworkException <: Exception 
    var::AbstractString 
end 

mutable struct ICException <: Exception 
    var::AbstractString 
end 

mutable struct ControlException <: Exception 
    var::AbstractString 
end 

mutable struct MissingDataException <: Exception 
    var::AbstractString 
end 

Base.showerror(io::IO, e::CFLException) = print(io, "CFL condition fails, time-step too large for component ", e.var, "!")

Base.showerror(io::IO, e::NetworkException) = print(io, "Network topology incompatible with simulator; ", e.var, "!")

Base.showerror(io::IO, e::ICException) = print(io, "Initial condition missing: ", e.var, "!")

Base.showerror(io::IO, e::ControlException) = print(io, "Control values: ", e.var, "!")

Base.showerror(io::IO, e::MissingDataException) = print(io, "Data: ", e.var, "missing!")

ref(ts::TransientSimulator) = ts.ref
ref(ts::TransientSimulator, key::Symbol) = ts.ref[key]
ref(ts::TransientSimulator, key::Symbol, id::Int64) = ts.ref[key][id]
ref(ts::TransientSimulator, key::Symbol, id::Int64, field) = ts.ref[key][id][field]

params(ts::TransientSimulator) = ts.params
params(ts::TransientSimulator, key::Symbol) = ts.params[key]

nominal_values(ts::TransientSimulator) = ts.nominal_values
nominal_values(ts::TransientSimulator, key::Symbol) = ts.nominal_values[key]

initial_pipe_mass_flow(ts::TransientSimulator, id::Int64) = 
    ts.initial_conditions[:pipe]["mass_flow"][id]

initial_pipe_pressure(ts::TransientSimulator, id::Int64) = 
    get(ts.initial_conditions[:pipe]["pressure"], id, false)

initial_nodal_pressure(ts::TransientSimulator, id::Int64) = 
    ts.initial_conditions[:node][id]

has_compressor_initial_flows(ts::TransientSimulator) = 
    !isempty(ts.initial_conditions[:compressor])

initial_compressor_flow(ts::TransientSimulator, id::Int64) = 
    ts.initial_conditions[:compressor][id]

function control(ts::TransientSimulator,
    key::Symbol, id::Int64, t::Real)::Tuple{CONTROL_TYPE,Float64}
    (key == :node) && (return get_nodal_control(ts, id, t))
    (key == :compressor) && (return get_compressor_control(ts, id, t))
    throw(ControlException("control available only for nodes and compressors"))
    return CONTROL_TYPE::unknown_control, 0.0
end

get_eos_coeffs(ts::TransientSimulator) = ts.pu_eos_coeffs(nominal_values(ts), params(ts))
get_pressure(ts::TransientSimulator, density) = ts.pu_density_to_pu_pressure(density, nominal_values(ts), params(ts))
get_density(ts::TransientSimulator, pressure) = ts.pu_pressure_to_pu_density(pressure, nominal_values(ts), params(ts))

TOL = 1.0e-5

function get_nodal_control(ts::TransientSimulator,
    id::Int64, t::Real)::Tuple{CONTROL_TYPE,Float64}
    if !haskey(ts.boundary_conditions[:node], id)
        return flow_control, 0.0
    end

    control_type = ts.boundary_conditions[:node][id]["control_type"]

    if control_type == flow_control
        val = ts.boundary_conditions[:node][id]["spl"](t - ts.params[:dt]/2)
    else
        val = ts.boundary_conditions[:node][id]["spl"](t)
    end

    return control_type, val
end

function get_compressor_control(ts::TransientSimulator,
    id::Int64, t::Real)::Tuple{CONTROL_TYPE,Float64}
    control_type, val = ts.boundary_conditions[:compressor][id]["spl"](t)

    if CONTROL_TYPE(control_type) == flow_control
        control_type_prev, val = ts.boundary_conditions[:compressor][id]["spl"](t - ts.params[:dt]/2)
        @assert control_type == control_type_prev
    end

    return CONTROL_TYPE(control_type), val
end

struct CompressorControl
    time::Vector{Float64}
    control_type::Vector{Any}
    value::Vector{Float64}
end

function evaluate(sp::CompressorControl, t::Real)::Tuple{Any,Float64}
    # @assert (t >= sp.time[1] - TOL)
    time = sp.time
    control_type = sp.control_type
    value = sp.value
    (isapprox(t, time[1], rtol=TOL)) || (t<= time[1]) && (return control_type[1], value[1])
    (isapprox(t, time[end], rtol=TOL)) || (t >= time[end]) && (return control_type[end], value[end])
    index = findlast(x -> t >= x, time)
    (isapprox(t, time[index], rtol=TOL)) && (return control_type[index], value[index])
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

# call synonym for CompressorControl evaluation
(spl::CompressorControl)(t::Real) = evaluate(spl, t)

@enum CONTROL_TYPE begin
    c_ratio_control = 0
    discharge_pressure_control = 1
    flow_control = 2
    pressure_control = 10
    unknown_control = 100
end

get_barglyphs() = 
    BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',)