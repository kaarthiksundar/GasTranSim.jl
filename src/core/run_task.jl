function _execute_task!(task_func::Function, ts::TransientSimulator, key_array::Vector, run_type::Symbol)
    if (run_type == :async)
        @sync for key in key_array
            @async task_func(ts, key)
        end
    else
        for key in key_array
            task_func(ts, key)
        end
    end
    return
end
