# GasTranSim.jl Documentation

```@meta
CurrentModule = GasTranSim
```

## Overview

Before illustrating the call  to invoke the simulator, we go over  the various input data files  and their contents with reference to the simple case of  a single pipe. These input files are contained in the directory `examples/data/model1pipe_fast_transients/`


## Description of network (network.json)

A single pipe gives rise to a network structure of two vertices (ends of the pipe) and one edge (pipe). This topological information is encoded in `network.json`. Information about compressors is included here if present and nodes which have pressure specified (slack-nodes) are indicated through a boolean flag.

	{
	    "nodes": {
	        "1": {
	            "node_id": 1,
	            "node_name": "n1",
	            "x_coord": -1,
	            "y_coord": 0,
	            "min_pressure": 3000000,
	            "max_pressure": 6000000,
	            "min_injection": 0,
	            "max_injection": 1000,
	            "slack_bool": 1
	        },
	        "2": {
	            "node_id": 2,
	            "node_name": "n2",
	            "x_coord": -0.0421,
	            "y_coord": 0,
	            "min_pressure": 3000000,
	            "max_pressure": 6000000,
	            "min_injection": -1000,
	            "max_injection": 1000,
	            "slack_bool": 0
	        }
	    },
	    "pipes": {
	        "1": {
	            "pipe_id": 1,
	            "pipe_name": "p1",
	            "from_node": 1,
	            "to_node": 2,
	            "diameter": 0.9144,
	            "length": 20000,
	            "friction_factor": 0.01,
	            "disc_seg": 0
	        }
	    }
	}	


Note that at quantities like diameter, length and pressure have numbers associated, but we do not know what  units they refer to. The required information to interpret these numbers is specified in another input file  which we now describe.

## Problem parameters (params.json)
This file records the parameters of the problem and the units in which they are specified.


	{
    "simulation_params": {
        "Temperature (K):": 239.11,
        "Gas specific gravity (G):": 0.6,
        "Specific heat capacity ratio": 1.4,
        "units (SI = 0, standard = 1)": 0.0,
        "Initial time": 0,
        "Final time": 3600,
        "Discretization time step": 1,
        "Courant number (must be between 0 and 1, recommended value is 0.9)": 0.90,
        "Output dt": 1.0,
        "Output dx": 1000.0,
        "Save final state": 1.0
    	}
	}


Thus, the fixed time step for the simulation, the time step to be used for  saving the results as output if indicated are all specified here.

## Boundary conditions (bc.json)

For the case of a single pipe, conditions have to be specified at the two ends. 
The `network.json` file already indicated  that the left end is a slack node, where pressure will be specified for all times.  Gas withdrawal is specified as a function of time at the other end.
The list of discrete time instants and the corresponding value are specified in the json file  and these values are interpolated when necessary.

	{
	    "boundary_nonslack_flow": {
	        "2": {
	            "time": [
	                0,
	                599,
	                600,
	                1799,
	                1800,
	                86400
	            ],
	            "value": [
	                0,
	                0,
	                787.63, 
	                787.63,
	                78.76,
	                78.76
	            ]
	        }
	    },
	    "boundary_pslack": {
	        "1": {
	            "time": [
	                0,
	                86400
	            ],
	            "value": [
	                6.5e6,
	                6.5e6
	            ]
	        }
	    }
	}



## Initial condition (ic.json)

The initial condition should specify the pressure and mass flux throughout the network at the initial time instant. For this case in particular, the initial nodal pressures and the mass flux along the length of the pipe should be specified.

	{

	    "initial_nodal_flow": {
	        "1": 0,
	        "2": 0
	    },
	    "initial_nodal_pressure": {
	        "1": 6.5e6,
	        "2": 6.5e6
	    },
	    "initial_pipe_flow": {
	        "1": 0
	    },
	    "initial_pipe_pressure_in": {
	        "1": 6.5e6
	    },
	    "initial_pipe_pressure_out": {
	        "1": 6.5e6
	    }
	}

Note that here the initial conditions correspond to a steady flow which is why the flow is the same throughout the pipe. In general, the initial condition will have flows that are spatially varying, in which case instead of specifying a single scalar, the mass flux at discrete lengths along the pipe will be recorded, analogous to what we saw above with `bc.json` for temporal variation. We shall deal with such a case subsequently.

