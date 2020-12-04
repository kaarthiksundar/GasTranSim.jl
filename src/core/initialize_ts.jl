function initialize_simulator(file::AbstractString)::TransientSimulator
    data = parse_json(file)
    return initialize_simulator(data)
end 

function initialize_simulator(data::Dict{String,Any})::TransientSimulator
    params, nominal_values = process_data!(data)
    make_per_unit!(data, params, nominal_values)
    ref = build_ref(data, ref_extensions= [add_pipe_info_at_nodes!, add_compressor_info_at_nodes!,
        add_current_time_to_ref!]) 
    ref[:current_time] = params[:t_0]
    bc = build_bc(data)
    pu_pressure_to_pu_density = x -> x 
    pu_density_to_pu_pressure = x -> x 
    pu_dp_drho = x-> x^0

    return TransientSimulator(data, 
        ref, 
        Dict{String,Any}(), 
        nominal_values, 
        params, 
        bc, 
        pu_pressure_to_pu_density, 
        pu_density_to_pu_pressure,
        pu_dp_drho
        ) 
end 



function add_grid_to_ref!(ts::TransientSimulator)
    for (key, pipe) in ts.ref[:pipe]
        n = ceil(Int64, ts.ref[:pipe][key]["length"] / (1 * ts.params[:dt]) ) 
        ts.ref[:pipe][key]["num_discretization_points"] = n
        ts.ref[:pipe][key]["dx"] = ts.ref[:pipe][key]["length"] / (n-1)
        ts.ref[:pipe][key]["density_profile"] = zeros(Float64, n)
        ts.ref[:pipe][key]["mass_flux_profile"] = zeros(Float64, n+1)
    end
end


function advance_current_time!(ts::TransientSimulator, t::Real)
    ts.ref[:current_time] += t
end

##pipes
function advance_mass_flux_internal!(ts::TransientSimulator)
    fun(a::Real, y::Real) = sign(y)  * ( -1 + sqrt( 1 + 4*a*abs(y) ) )/ (2*a)

    for (key, pipe) in ts.ref[:pipe]
        rho = ts.ref[:pipe][key]["density_profile"]
        phi = ts.ref[:pipe][key]["mass_flux_profile"]
        n = ts.ref[:pipe][key]["num_discretization_points"]

        for i = 2:n
            a = ts.params[:dt] * ts.ref[:pipe][key]["friction_factor"] / (rho[i] + rho[i-1])
            y = phi[i] - ( ts.params[:dt] / ts.ref[:pipe][key]["dx"] ) * (get_pressure(ts, rho[i]) - get_pressure(ts, rho[i-1])) - a * phi[i] * abs(phi[i])
            phi[i] = fun(a, y)
        end
        
        # update field. Can also loop over all pipes again instead using function below
        ts.ref[:pipe][key]["from_minus_mass_flux"] = phi[2]
        ts.ref[:pipe][key]["to_minus_mass_flux"] = phi[n]
    end

end

function update_field_mass_flux_minus!(ts::TransientSimulator)
    for (key, pipe) in ts.ref[:pipe]
        ts.ref[:pipe][key]["from_minus_mass_flux"] = ts.ref[:pipe][key]["density_profile"][2]
        ts.ref[:pipe][key]["to_minus_mass_flux"] = ts.ref[:pipe][key]["density_profile"][ts.ref[:pipe][key]["num_discretization_points"]]        
    end
end

function advance_density_internal!(ts::TransientSimulator)

    for (key, pipe) in ts.ref[:pipe]
        rho = ts.ref[:pipe][key]["density_profile"]
        phi = ts.ref[:pipe][key]["mass_flux_profile"]
        n = ts.ref[:pipe][key]["num_discretization_points"]

        for i = 2:n-1
        rho[i] += ( ts.params[:dt] / ts.ref[:pipe][key]["dx"] ) *  (phi[i] - phi[i+1])
        end

    end

end



## Nodes
function _set_pressure_at_vertex!(vertex_key::Int64, val::Float64, ts::TransientSimulator)
     if ts.ref[:nodes][vertex_key]["is_updated"] == False
        ts.ref[:nodes][vertex_key]["pressure"] = val
        ts.ref[:nodes][vertex_key]["is_updated"] = True
    end
end

function _set_pressure_at_vertex_across_compressors!(vertex_key::Int64, prs::Float64, ts::TransientSimulator)
    
    if !isempty(ts.ref[:incoming_compressors][vertex_key])
        for ci in ts.ref[:incoming_compressors][vertex_key]
            ctrl_type, val = control(ts, :compressor, ci, ts.ref[:current_time])           
            i = ts.ref[:compressors][ci]["fr_node"]

            if ctrl_type == c_ratio_control              
                _set_pressure_at_vertex!(i, prs/val, ts)
                # do an update for all compressors together instead of now
                # ts.ref[:compressors][ci]["discharge_pressure"] = prs
            else if ctrl_type == discharge_pressure_control
                @error "compressor discharge pressure control coincides with nodal pressure control"
            end
        end
    end

    if !isempty(ts.ref[:outgoing_compressors][vertex_key])
        for ci in ts.ref[:outgoing_compressors][vertex_key]
            ctrl_type, val = control(ts, :compressor, ci, ts.ref[:current_time])       
            i = ts.ref[:compressors][ci]["to_node"]

            if ctrl_type == c_ratio_control            
                _set_pressure_at_vertex!(i, prs * val, ts)
                # ts.ref[:compressors][ci]["discharge_pressure"] = prs * val
            else if ctrl_type == discharge_pressure_control
                _set_pressure_at_vertex!(i, val, ts)
                # ts.ref[:compressors][ci]["c_ratio"] = val/prs
            end
        end
    end

end


function _solve_for_pressure_at_vertex!(vertex_key::Int64, val::Float64, ts::TransientSimulator)



end

function advance_pressure_vertex!(ts::TransientSimulator)

    for (key, junction) in ts.ref[:node]
        
        if ts.ref[:nodes][key]["is_updated"] == True
            continue
        end
        
        # p(t), but q(t - dt/2)
        ctrl_type, val = control(ts, :node, key, ts.ref[:current_time])
        
        if ctrl_type == pressure_control
            _set_pressure_at_vertex!(key, val, ts)
            _set_pressure_at_vertex_across_compressors!(key, val, ts)
            
        else if ctrl_type == flow_control
            _solve_for_pressure_at_vertex!(key, val, ts)
            println(ctrl_type, typeof(ctrl_type))
        end
    end

# else 
 #loop over list of incoming pipes, outgoing pipes
 # loop over list of incoming compressors
   # for each compressor, take action dependent on type of compressor
 # similarly loop over list of outgoing compressors
 
end

