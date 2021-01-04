function _execute_task!(task_func::Function, ts::TransientSimulator, key_array::Vector, run_type::Symbol)
    num_jobs = length(key_array)
    jobs = Channel{Any}(num_jobs)
    @async begin
        for key in key_array
            put!(jobs, key)
        end
        close(jobs)
    end

    if (run_type == :async)
        @sync for key in jobs
            @async task_func(ts, key)
        end
    else
        for key in jobs
            task_func(ts, key)
        end
    end
    return
end
