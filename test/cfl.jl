using GasTranSim

function main()
    folder = "./test/data/8-node-steady/"
    dt_factor = 1.0 

    while true
        try 
            data = get_data(folder) 
            delta_t = data["simulation_params"]["Discretization time step"]
            data["simulation_params"]["Discretization time step"] = delta_t * dt_factor

            # ideal run 
            ts = initialize_simulator(data; eos=:ideal)
            run_simulator!(ts)    
            break
            
        catch err 
                
            if isa(err, CFLException)
                dt_factor *= 0.5
                println("CFL condition failed. Reducing Δt by half")
            end 
            
        end
    end 
end 

main()





