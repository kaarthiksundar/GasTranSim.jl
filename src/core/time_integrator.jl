function set_up_grid_from_CFL!(dt::Real, ts::TransientSimulator, integrator::ExplicitStaggeredSpaceTimeIntegrator)

	for (key, pipe) in ts.ref[:pipe]
		n = ceil( pipe["length"] / 1 * dt)
		integrator.pipe_edge[key][:numpts] = n
		integrator.pipe_edge[key][:dx] = pipe["length"]/(n-1)
		integrator.pipe_edge[key][:rho] = zeros(Float64, n)
		integrator.pipe_edge[key][:phi] = zeros(Float64, n+1)

	end
end


function advance_phi_internal!(dt::Real, ts::TransientSimulator, integrator::ExplicitStaggeredSpaceTimeIntegrator)

	for pipe in integrator.pipe_edge
		for i = 2: pipe[:numpts]
			pipe[:phi][i] +=
		end
	end

end

function advance_rho_internal!(dt::Real, ts::TransientSimulator, integrator::ExplicitStaggeredSpaceTimeIntegrator)

	for pipe in integrator.pipe_edge
		for i = 2: pipe[:numpts] -1
			pipe[:rho][i] += (dt/pipe[:dx]) *  (pipe[:phi][i] - pipe[:phi][i+1])
		end
	end

end

function advance_rho_vertex(dt)