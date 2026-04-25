function initialize_pipe_grid!(ts::TransientSimulator, method::Symbol)
    return initialize_pipe_grid!(ts, Val(method))
end

function initialize_pipe_grid!(::TransientSimulator, ::Val{method}) where {method}
    throw(
        ArgumentError(
            "Unknown spatial discretization method :$method. Supported methods: $(_supported_methods_str())",
        ),
    )
end

function initialize_pipe_state!(ts::TransientSimulator, method::Symbol)
    return initialize_pipe_state!(ts, Val(method))
end

function initialize_pipe_state!(::TransientSimulator, ::Val{method}) where {method}
    throw(
        ArgumentError(
            "Unknown pipe-state initialization method :$method. Supported methods: $(_supported_methods_str())",
        ),
    )
end
