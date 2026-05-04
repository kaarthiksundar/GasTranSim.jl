const _SUPPORTED_METHODS = (:explicit_staggered_grid, :explicit_staggered_grid_new, :implicit_parabolic, :explicit_hyperbolic)

supported_methods() = collect(_SUPPORTED_METHODS)

function _supported_methods_str()
    return join(string.(supported_methods()), ", ")
end

function assert_supported_method!(method::Symbol)
    method in _SUPPORTED_METHODS && return
    throw(
        ArgumentError(
            "Unknown method :$method. Supported methods: $(_supported_methods_str())",
        ),
    )
end

function validate_method_contract!(method::Symbol)
    assert_supported_method!(method)

    missing = String[]
    hasmethod(initialize_pipe_grid!, Tuple{TransientSimulator,Val{method}}) ||
        push!(missing, "initialize_pipe_grid!(ts, ::Val{$method})")
    hasmethod(initialize_pipe_state!, Tuple{TransientSimulator,Val{method}}) ||
        push!(missing, "initialize_pipe_state!(ts, ::Val{$method})")
    hasmethod(step!, Tuple{TransientSimulator,Val{method},Symbol}) ||
        push!(missing, "step!(ts, ::Val{$method}, run_type)")
    

    isempty(missing) && return
    throw(
        ArgumentError(
            "Method :$method is missing required implementations: $(join(missing, "; "))",
        ),
    )
end
