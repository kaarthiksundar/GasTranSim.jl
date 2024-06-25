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

## Nondimensionalized Equations

Let us choose the nominal quantities $p_0, \rho_0, v_0, \phi_0, t_0, l_0,f_0, A_0$ for pressure, density, velocity, mass flux, time, length, mass flow and cross-sectional area respectively.

Then we can rewrite the governing equations in nondimensional variables $\bar{\rho}, \bar{t},\bar{\phi}, \bar{x}, \dotsc$  through the substitutions $\rho = \rho_0 \bar{\rho}, t = t_0\bar{t}, \dotsc$.

Then 
$$\dfrac{\partial \bar{\rho}}{\partial \bar{t}} + \left( \dfrac{\phi_0t_0}{\rho_0l_0} \right)\dfrac{\partial \bar{\phi}}{\partial \bar{x}} = 0$$
and
$$\dfrac{\partial \bar{\phi}}{\partial \bar{t}} + \left( \dfrac{p_0t_0}{l_0\phi_0} \right)\dfrac{\partial \bar{p}}{\partial \bar{x}} = -\dfrac{\lambda}{2D}\left( \dfrac{\phi_0t_0}{\rho_0} \right) \dfrac{\bar{\phi}|\bar{\phi}|}{\bar{\rho}}$$

We use a *fixed* value for the speed of sound $a$, choose values for $l_0, p_0, \rho_0, v_0$ based on the data, and set $A_0 =1, \phi_0 = \rho_0 v_0, t_0 = l_0/v_0, f_0 = \phi_0 A_0$.

This choice of values reduces the governing equations to 
$$\dfrac{\partial \bar{\rho}}{\partial \bar{t}} + \dfrac{\partial \bar{\phi}}{\partial \bar{x}} = 0$$
$$\dfrac{\partial \bar{\phi}}{\partial \bar{t}} + \dfrac{\mathcal{C}}{\mathcal{M}^2}\dfrac{\partial \bar{p}}{\partial \bar{x}} = -\left( \dfrac{\lambda}{2D/l_0} \right) \dfrac{\bar{\phi}|\bar{\phi}|}{\bar{\rho}}$$

where $\mathcal{C} =  \dfrac{p_0}{\rho_0 a^2}$ and $\mathcal{M} = \dfrac{v_0}{a}$ are the Euler number and Mach number respectively.

The ideal gas equation, $p = \rho R_g T$ transforms to $\bar{\rho} = \mathcal{C} \bar{p}$ since $a = \sqrt{R_gT}$. In non-isothermal problems, we can define a nominal temperature $T_0$  and set $a=\sqrt{R_gT_0}$ in that case to get $\bar{\rho} = \dfrac{\mathcal{C}T_0}{T} \bar{p}$.

For a non-ideal gas, assuming the CNGA equation of state we have $p\cdot(b_1 + b_2 p) = \rho R_g T$ which simplifies to the expression $\bar{\rho} = \left ( \bar{b_1} + \bar{b_2}\bar{p} \right ) \bar{p}$, where $\bar{b_1} = \mathcal{C}b_1, \bar{b_2} = \mathcal{C}p_0 b_2$. 

The speed of sound $a$ in a gas  is given by $a^2 = \dfrac{d p}{d \rho}$. 
For an ideal gas, we get $a^2 = R_g T$, or in nondimensional terms $\dfrac{d\bar{p}}{d\bar{\rho}} = 1/\mathcal{C} = \bar{a}^2$.

For a simple CNGA equation, $p\cdot(b_1 + b_2 p) = \rho R_g T$, we can calculate the sound speed from $\dfrac{d \rho}{d p}$ since $\dfrac{d p}{d \rho} \cdot \dfrac{d \rho}{d p} = 1$. Thus  $a_{CNGA}^2 = \dfrac{R_g T}{b_1 + 2b_2 p}$. Since $b_1 > 1, b_2 > 0, p > 0$, $a_{CNGA} < a$.
We can reach the same conclusion from the fact that  $\dfrac{d \bar{p}}{d \bar{\rho}} = \dfrac{1}{\bar{b_1} + 2\bar{b_2}\bar{p}}$ and $\bar{b_1} > \mathcal{C}, \bar{b_2} > 0$.

Let us record the CFL condition and its consequence next. Usually, it would be stated as $a \Delta {t}/ \Delta {x} \leq 1$, but to be safe, we consider $ a \Delta {t}/ \Delta {x} \leq k$ for $k= 0.9$. 
The dimensionless form of this condition is $\Delta \bar{t}/ (\mathcal{M} \Delta \bar{x}) \leq k$.


The CFL condition  was derived in the paper assuming the ideal gas equation but we use it for non-ideal gas as well since we have relaxed $k=0.9$.  For the particular choice $v_0 = a$, we  get $\mathcal{C} = \mathcal{M}  = 1$.

Setting $\Delta \bar{x} = \bar{L}/m$,  we can find the number of pipe segments $m$.


