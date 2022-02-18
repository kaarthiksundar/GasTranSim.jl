# GasTranSim.jl Documentation

```@meta
CurrentModule = GasTranSim
```

## Governing Equations and Discretization
The governing equations and the equation of state are formulated in terms of the density $\rho$, pressure $p$, and mass-flux $\phi$ for a pipe with given diameter $D$ and known friction factor $\lambda$ as

$$\partial_t \rho + \partial_x \phi = 0, \partial_t \phi + \partial_x p = -\frac{\lambda}{2D}\frac{\phi |\phi|}{\rho}, p = Z(p, T)RT\rho$$

This system is augmented with given initial conditions and boundary conditions. For a network, initial conditions specify $\rho, p, \phi$ everywhere at initial time $t = t_0$.
For  times $t > t_0$, at every junction in the network, the withdrawal (mass-flow per unit time) is assumed to be zero, unless specified otherwise, while the pressure/density is specified at \emph{at least} one junction. These are the boundary conditions for the problem.

For each compressor, its action is assumed to be given through specification of one of the following quantities for all time - (i) the compressor ratio (ii) the compressor delivery pressure (iii) the mass-flow through the compressor. 

In order to solve the problem, we discretize the governing equations using the explicit, second order, staggered finite difference  discretization scheme proposed by
* V. Gyrya and A. Zlotnik (2019). An explicit staggered-grid method for numerical simulation of large-scale natural gas pipeline networks. ([doi:10.1016/j.apm.2018.07.051](https://doi.org/10.1016/j.apm.2018.07.051))



