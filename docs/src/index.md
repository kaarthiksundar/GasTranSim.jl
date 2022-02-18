# GasTranSim.jl Documentation

```@meta
CurrentModule = GasTranSim
```

## Overview

GasTranSim.jl is a Julia/JuMP package for simulating transient flow of a gas through a given network of gas pipelines which may include compressors/boosting stations with defined operation. 
The simulator allows us to predict the pressure and mass-flux (mass-flow per unit time for unit cross-sectional area) throughout the pipeline network for a given  time interval.
Moreover, it allows the user to save the state of the simulation at a particular time instant, and restart the simulation seamlessly from the saved state.

## Description of network

The gas pipeline network may be thought of as a graph, where the edges are pipelines/compressors, and vertices are junctions of different pipes and/or compressors.
The pipeline junctions in the network can have specified gas supply/withdrawal (mass flow per unit time) or a specified pressure. The operation of a compressor element is defined by specifying one of the following: (i) compressor ratio (r > 1, ratio of pressures across the ends of the compressor) (ii) delivery pressure (iii) mass flow per unit time through the compressor.

## Installation Guide

To use GasTranSim, first [download and install](https://julialang.org/downloads/) Julia or open up a remote notebook at [JuliaBox](https://www.juliabox.com/) or similar services.

This version of GasTranSim is compatible with Julia 1.0 and later.

From Julia REPL, GasTranSim is installed by using the built-in package manager:
```julia
import Pkg
Pkg.add("GasTranSim")
```

## Unit Tests
To run the tests in the package, run the following command within the Julia REPL after installing the package.

```julia
import Pkg
Pkg.test("GasTranSim")
```