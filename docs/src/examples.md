# GasTranSim.jl Documentation

```@meta
CurrentModule = GasTranSim
```

## Overview

Some example problems have been solved to illustrate the use of this simulator.  Two of the problems are based on the case of a single pipe, while the other two  are for  networks with 8 and 30 nodes respectively.

## Single pipe and 8 node network example

These three examples are taken from 
* V. Gyrya and A. Zlotnik (2019). An explicit staggered-grid method for numerical simulation of large-scale natural gas pipeline networks. ([doi:10.1016/j.apm.2018.07.051](https://doi.org/10.1016/j.apm.2018.07.051))


- The directory `1-pipe-fast-transients` contains the input data and source files related to Section 5.2 in the article.

- The directory `1-pipe-slow-transients` pertains to Section 5.3 in the paper.

- The directory `8-node` replicates some results of the problem discussed in Section 6 of the paper.


## 30 node network example

The directory  `model30` has the input data and source files to demonstrate that if we start with a steady state initial condition,  the same steady state solution can be recovered from the transient simulation.

## Running the examples

Each of the four examples has its own file that can be invoked from the root directory as follows 
```julia
import Pkg 
Pkg.activate("examples/")
include("examples/1-pipe-fast.jl")
```







