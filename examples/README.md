# GasTranSim.jl examples

This folder contains examples. There are three examples: 

1. ``1-pipe-fast-runs.jl``
2. ``1-pipe-slow-runs.jl`` 
3. ``8-node=runs.jl`` 

All three examples are taken from the paper "An explicit staggered-grid method for numerical simulation of large-scale natural gas pipeline networks" by Vitaliy Gyrya, Anatoly Zlotnik. The arXiv version of the paper can be found in the [arXiv link](https://arxiv.org/abs/1803.00418). The results of the runs reproduce the results in the paper (for more details see the examples page in the Documentation)

To run the three examples use the following command: 

1. ``julia --project=examples/ 1-pipe-fast-runs.jl`` 
2. ``julia --project=examples/ 1-pipe-slow-runs.jl``
3. ``julia --project=examples/ 8-node-runs.jl``

The output figures are stored in the ``examples/tmp/`` directory.
