function step!(ts::TransientSimulator, method::Symbol, run_type::Symbol)
    return step!(ts, Val(method), run_type)
end

function step!(ts::TransientSimulator, ::Val{:explicit_staggered_grid}, run_type::Symbol)
    explicit_staggered_grid_step!(ts, run_type)
    return
end

function step!(ts::TransientSimulator, ::Val{:implicit_parabolic}, run_type::Symbol)
    implicit_parabolic_step!(ts, run_type)
    return
end

function step!(::TransientSimulator, ::Val{method}, ::Symbol) where {method}
    throw(
        ArgumentError(
            "Unknown time-integration method :$method. Supported methods: $(_supported_methods_str())",
        ),
    )
end
