# GasTranSim.jl Documentation

```@meta
CurrentModule = GasTranSim
```

## Overview

GasTranSim.jl is a Julia/JuMP package for 

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