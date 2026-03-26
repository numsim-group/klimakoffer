+++
title = "Milestone 6"
hascode = true
rss = "Description"
rss_title = "Milestone 6"
rss_pubdate = Date(2022, 5, 1)

tags = ["discrete ebm"]
+++

# Milestone 6 - Fully Discrete EBM

\toc

## Recap Semi-discretization

We derived the final form of the $2D$ EBM with heat diffusion in [milestone 5 in the section about planetary heat transfer](/milestone5/milestone5_heat_transfer/#choice_of_the_diffusion_coefficient_in_the_2d_ebm)
$$
C(x) \partialderiv{T}{t} + A(CO_2) + B T - \Nabla \cdot (D\Nabla T) = S_{sol}(x,t),
$$
where we here use the generic version in Cartesian coordinates for simplicity. 

The so-called semi-discretization of the PDE (discretization in spatial coordinates), derived in [milestone 5](/milestone5/milestone5_spatial_discretization/#semi-discrete_ebm), is given for a grid node $(j,i)$ as
$$
\deriv{T_{j,i}}{t} = 
\underbrace{
    \frac{L_{j,i}}{C_{j,i}} - \frac{B}{C_{j,i}} T_{j,i}
 }_{R_{j,i}(T)}
+
\underbrace{
    \frac{1}{C_{j,i}} \left(S_{sol, j, i}(t) - A\right)
}_{F_{j,i}},
$$
where $L_{j,i}$ denotes the spatial discretization of the diffusion operator $\Nabla \cdot (D\Nabla T)$. Recall, that for the spatial discretization, due to the pole problem of the spherical coordinates, a special treatment of the diffusion operator in the pole region is necessary. This form is called semi-discretization, as only the spatial directions are discretized, i.e., the formulation does not depend on the spatial coordinates $x$ anymore, but still on the temporal coordinate $t$. 

The resulting problem is a (system) of ODEs that need to be integrated in time. In particular, we start with the equivalent [form that we derived](/milestone5/milestone5_stability_jacobians/#eqsystem_linear) taking into account the linearity of the diffusion operator and the linearity of the radiative heat transfer,
$$\label{eq:system_linear}
\dot{\mathbf{T}} = 
\underbrace{
    \mat{A} \mathbf{T}
}_{\mathbf{R}(\mathbf{T})}
+
\mathbf{F}(t).
$$
Due to the linearity of the semi-discretization in the unknown variable (temperature $\mathbf{T}$ in case of the EBM), we can represent the spatial EBM operator $\mathbf{R}(\mathbf{T})$ as a constant matrix times the vector of unknown temperatures, $\mathbf{R}(\mathbf{T}) = \mat{A} \mathbf{T}$. The matrix $\mat{A}\in \R^{\ndof\times\ndof}$ is the so-called Jacobian matrix of the spatial operator and its computation is described in [milestone 5](/milestone5/milestone5_stability_jacobians/#computing_jacobian_matrices_numerical_jacobian).

## Time Integration

We can now apply the two time-integration methods that we learned in [milestone 3](/milestone3/milestone3_odesolvers/#ode_solvers) to all equations of the linear system \eqref{eq:system_linear}.
For the evolution from time step $n$ ($t = t^n$) to time step $n+1$ ($t = t^{n+1}$), we discussed two options:

* If we use the forward Euler method, the fully discrete scheme reads as
    $$\label{eq:forward_euler}
    \frac{\mathbf{T}^{n+1} - \mathbf{T}^n}{\Delta t} = 
        \mat{A} \mathbf{T}^n
    +
    \mathbf{F}^n,
    $$
    and the solution at the time step $n+1$ can be obtained explicitly as
$$\label{eq:forward_euler2}
\mathbf{T}^{n+1} = \mathbf{T}^n + \Delta t \left(
\mat{A} \mathbf{T}^n
+
\mathbf{F}^n
\right),
$$
or eqvivalently
$$\label{eq:forward_euler3}
\mathbf{T}^{n+1} = \mathbf{T}^n + \Delta t \left(\mathbf{R}(\mathbf{T}^n)
+
\mathbf{F}^n
\right),
$$
i.e. in the case of the explicit method, no Jacobian computation is necessary, instead the operator can be used directly. Thus, per time step, the update with the explicit Euler method is very quick - but the extreme time step restriction for (von Neumann) stability of the fully-discrete scheme of the EBM is to restrictive in practise: although one explicit time step is cheap, the amount of time steps to reach the equilibrium is so high (due to the restrictive time step stability) that the overall simulation takes a very long time. 

* As a more efficient alternative, we discussed the backward Euler method, where the fully discrete scheme reads as
$$\label{eq:backward_euler}
\frac{\mathbf{T}^{n+1} - \mathbf{T}^n}{\Delta t} = 
    \mat{A} \mathbf{T}^{n+1}
+
\mathbf{F}^{n+1},
$$
and the solution at the time step $n+1$, $\mathbf{T}^{n+1}\in \R^{\ndof}$, must be obtained by solving the system of linear equations
$$\label{eq:backward_euler2}
\underbrace{
\left( \mat{I} - \Delta t\,\mat{A} \right) 
}_{\mat{M}}
\mathbf{T}^{n+1} = 
\underbrace{\mathbf{T}^n
+
\Delta t\,\mathbf{F}^{n+1}}_{\mathbf{Y}^{n+1}},
$$
where $\mat{I}\in \R^{\ndof\times\ndof}$ is the identity matrix, $\mat{M}\in \R^{\ndof\times\ndof}$ is the system matrix of the linear equation system, and $\mathbf{Y}^{n+1}\in \R^{\ndof}$ is the right-hand side vector of the linear system.

@@colbox-blue
**Remark:** We note that for the presented 2D EBM, the system matrix $\mat{M}$ is constant in time if $\Delta t$ is fixed. Hence, it can be computed once at the beginning of the simulation and does not change throughout the whole computation. 
This of course would change, if more complex processes are involved in the EBM.
For instance, feedback between geography and the resulting temperature field, or in the case of a non-linear radiation model (e.g. Stefan-Boltzmann law). The right-hand side vector $\mathbf{Y}^{n+1}$, on the other hand, changes for each time step and needs to be re-computed accordingly. 
@@

## Equilibrium simulation

Similar to the task in the pure ODE case, see [milestone 3](/milestone3/milestone3_odesolvers/#temporal_equilibrium_simulation), we are seeking to simulate the equilibrium of the temperature field throughout the year. 
We emphasize again that we are not looking for a single steady-state temperature field - due to the annual dependence of the solar forcing, no global steady state solution can be reached. 
Instead, it is possible to reach an annual equilibrium state, where the solution changes in time, but repeats itself every year, such that we call it an equilibrium. 

Again, the task is to define a proper stopping criterium when performing the temporal EBM simulations. One could again compute an annual average temperature and stop when the change of said temperature is smaller than a given threshold. One could also compare directly the full 2D temperature field throughout the year and define a proper norm to define a stopping criterium. 

In summary, analogous to the pure ODE case, we have to compare the yearly temperature solutions with each other, until the difference from one year to the next is smaller than a given tolerance. Among the different choices of norms to compute the yearly solutions, the simplest option is to compute the Euclidean norm of the data vectors in time.