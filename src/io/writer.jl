function write_output(ts::TransientSimulator;
    output_path::AbstractString="./", 
    output_file::AbstractString="output.json", 
    final_state_file::AbstractString="ic_final.json")

    solution = ts.sol 
    output_string = output_path * output_file
    output_dict = Dict(
        "time_step"  => solution["time_step"],
        "time_points" => solution["time_points"],
        "nodes" => solution["nodes"],
        "pipes" => solution["pipes"],
        "compressors" => solution["compressors"]
    )

    open(output_string, "w") do f 
        JSON.print(f, output_dict, 2)
    end 

    save_final_state = params(ts, :save_final_state)
    if save_final_state == 1
        final_state_string = output_path * final_state_file
        final_state_dict = Dict(
            key => value for (key, value) in solution["final_state"]
        )

        open(final_state_string, "w") do f 
            JSON.print(f, final_state_dict, 2)
        end
    end 
end 

function write_final_state(ts::TransientSimulator; 
    output_path::AbstractString="./", 
    final_state_file::AbstractString="ic_final.json")

    save_final_state = params(ts, :save_final_state) 
    (save_final_state != 1) && (return)
    solution = ts.sol
    final_state_string = output_path * final_state_file
    final_state_dict = Dict(
        key => value for (key, value) in solution["final_state"]
    )

    open(final_state_string, "w") do f 
        JSON.print(f, final_state_dict, 2)
    end
end 