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

## Time Integration

We start from the semi-discrete EBM that we derived in milestone 5.
In particular, we use the [form that we derived](/milestone5/milestone5_stability_jacobians/#eqsystem_linear) taking into account the linearity of the diffusion operator and the radiative heat transfer,
$$\label{eq:system_linear}
\dot{\mathbf{T}} = 
\underbrace{
    \mat{A} \mathbf{T}
}_{\mathbf{R}(\mathbf{T})}
+
\mathbf{F}(t).
$$

We can now apply the two time-integration methods that we learned in [milestone 3](/milestone3/milestone3_odesolvers/#ode_solvers) to all equations of the linear system \eqref{eq:system_linear}.
For the evolution from time step $n$ ($t = t^n$) to time step $n+1$ ($t = t^{n+1}$), we obtain two possibilities:

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
\right) ,
$$

* If we use the backward Euler method, the fully discrete scheme reads as
$$\label{eq:backward_euler}
\frac{\mathbf{T}^{n+1} - \mathbf{T}^n}{\Delta t} = 
    \mat{A} \mathbf{T}^{n+1}
+
\mathbf{F}^{n+1},
$$
and the solution at the time step $n+1$ must be obtained by solving the system of linear equations
$$\label{eq:backward_euler2}
\underbrace{
\left( \frac{1}{\Delta t} \mat{I} - \mat{A} \right) 
}_{\mat{M}}
\mathbf{T}^{n+1} = 
\frac{1}{\Delta t} \mathbf{T}^n
+
\mathbf{F}^{n+1},
$$
where $\mat{I}$ is the identity matrix and we introduce the system matrix $\mat{M}$.
