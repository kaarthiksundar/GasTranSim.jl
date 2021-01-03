using JSON 
using Dierckx

include("io/json.jl")
include("io/data_utils.jl")

include("unit_conversion/unit_convertor_utils.jl")
include("unit_conversion/to_si.jl")
include("unit_conversion/to_english.jl")
include("unit_conversion/to_pu.jl")
include("unit_conversion/unit_convertors.jl")

include("core/types.jl")
include("core/ref.jl")
include("core/bc.jl")
include("core/initialize_ts.jl")
include("core/time_integration.jl")
include("io/output.jl")



file = "./data/model8ts_3d.json";

ts = initialize_simulator(file)
add_grid_to_ref!(ts)
initialize_vertex!(ts)
initialize_pipes!(ts)
out_int = initialize_output_struc(ts)

dt = ts.params[:dt]


run_type = :async
while ts.ref[:current_time] < 1e4*dt # ts.params[:t_f]

	
	advance_current_time!(ts, dt)


	#  if current_time is where some disruption occurs, modify ts.ref now
    
	advance_density_internal!(ts, run_type) #(n+1) level
	
	advance_pressure_mass_flux_vertex!(ts, run_type) #pressure (n+1), flux (n+1/2)

	advance_mass_flux_internal!(ts, run_type) # (n+1 + 1/2) level

	#  if current_time is one where output needs to be saved, check and do now
	update_output_struc!(ts, out_int)

	# err_arr = []
	# for i = 1:8
		
	# 	push!(err_arr, abs(ts.ref[:node][i]["pressure_previous"] - ts.ref[:node][i]["pressure"]) )
		
	# end
	# @show err_arr


end


# flux profile, density profiles saved to restart time marching
#out = create_output(ts, out_int)






 