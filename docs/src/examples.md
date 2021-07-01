# GasTranSim.jl Documentation

```@meta
CurrentModule = GasTranSim
```

## Overview

Some example problems have been solved to illustrate the use of this simulator.  Two of the problems are based on the case of a single pipe, while the other two  are for  networks with 8 and 30 nodes respectively.

## Single pipe and 8 node network example

These three examples are taken from 
* V. Gyrya and A. Zlotnik (2019). An explicit staggered-grid method for numerical simulation of large-scale natural gas pipeline networks. ([doi:10.1016/j.apm.2018.07.051](https://doi.org/10.1016/j.apm.2018.07.051))


- The directory `model1pipe_fast_transients` contains the input data and source files related to Section 5.2 in the article.

- The directory `model1pipe_slow_transients` pertains to Section 5.3 in the paper.

- The directory `model8_paper_VG_AZ` replicates some results of the problem discussed in Section 6 of the paper.


## 30 node network example

The directory  `model30` has the input data and source files to demonstrate that if we start with a steady state initial condition,  the same steady state solution can be recovered from the transient simulation.







