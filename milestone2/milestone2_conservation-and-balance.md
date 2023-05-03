+++
title = "Milestone 2"
hascode = true
rss = "Description"
rss_title = "Milestone 2"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 2 - Conservation and balance
\toc

## Conservation laws

Many processes in nature and engineering can be modeled (described) with a simple principle: the principle of **conservation**. Consider a closed domain $\Omega\in\mathbb{R}^n$ and a quantity $u=u(x,t)\in\mathbb{R}$ that is defined for all $x\in\Omega$ and $t\geq 0$. The function $u(x,t)$
typically describes a physical quantity such as mass, momentum, or energy.

We are interested in modeling the temporal evolution of $u$ (change in time $t$). Observations of nature of the behaviour of such quantities lead to the
following simple principle:
@@colbox-blue
**Observation:** The temporal change of $u(x,t)$ in a sub-domain $\omega\subset\Omega$ is equal to the amount
that gets generated or destroyed inside of $\omega$ in addition to the flux balance
inwards or outwards through the surface/boundary $\partial\omega$.
@@

The quantity of $u$ changes, if: 

(i) There is a positive/negative source inside of the sub-domain $\omega$

(ii) There is a positive/negative flux balance through the boundary of the sub-domain $\partial\omega$

@@colbox-blue
**Famous example:** Bathtub
\fig{/assets/milestone2/Bathtub.png}

The amount of water in the bathtub (\col{blue}{blue}) changes due to the sources (\col{green}{green}) and sinks (\col{red}{red}).
@@

## Formulation of the mathematical problem

Consider a sub-domain $\omega$ with outward pointing normal vector $\vec{n}$:

\fig{/assets/milestone2/domain.png}

$$
\underbrace{\partialderiv{}{t} \int_{\omega} u \d x}_{\substack{\text{temporal\,change} \\ \text{of}\,u\,\text{in}\,\omega}} = 
- \underbrace{\oint_{\partial \omega} \vec{f} \cdot \vec{n} \d s}_{\substack{\text{change~through} \\ \text{surface~flux}\,\vec{f}}} 
+ \underbrace{\int_{\omega} S(u,x,t) \d x}_{\substack{\text{change~through} \\ \text{source}~S(u,x,t)}} 
$$

@@colbox-blue
**Remark 1:** We call this the integral formulation of the problem.
@@

@@colbox-blue
**Remark 2:** The flux function $\vec{f}$ typically depends on the solution $u$ and the variables $x,t$, i.e., $\vec{f} = \vec{f}(u,x,t)$.
@@

If we make the mathematical assumption that the function $u(x,t)$ is sufficiently smooth (such that we can take the derivatives in space
and time), we can apply the Gauss integral theorem to the surface integral
\begin{align}
{\int_{\omega} \partialderiv{u}{t} \d x}
+ \underbrace{\int_{\omega} \Nabla \cdot \vec{f} \d x}_{=\int_{\omega} \partialderiv{f_1}{x_1} + \partialderiv{f_2}{x_2} + \cdots + \partialderiv{f_n}{x_n} \d x}
= \int_{\omega} S(u,x,t) \d x
\end{align}

As we have chosen an arbitrary sub-domain $\omega$ with no special properties, the equation needs to hold for all choices of $\omega\subset\Omega$. This
can only be true if the integrands balance out to zero, hence
$$
\frac{\partial u}{\partial t} + \vec{\nabla}\cdot\vec{f} = S,\quad x\in\Omega,\,\,t\in\mathbb{R}^+.
$$

@@colbox-blue
**Remark 3:** This is a partial differential equation (PDE) in space and time and is typically called a balance law.
@@

@@colbox-blue
**Remark 4:** For many processes, the source term can be neglected, i.e., we can choose $S(u,x,t) = 0$. The resulting PDE describes processes, where the quantity $u$ is neither destroyed nor generated, but is only changed by fluxes. Such PDEs are called conservation laws
$$
\frac{\partial u}{\partial t} + \vec{\nabla}\cdot\vec{f} = 0,\quad x\in\Omega,\,\,t\in\mathbb{R}^+.
$$
@@

Depending on the process we want to model, we need to chose the quantity of interest $u(x,t)$, a model for the flux $\vec{f}(u,x,t)$  and a model for the source term $S(u,x,t)$. There are 
many examples of such models, e.g., mass and momentum conservation in fluid mechanics. We consider here as an example the so-called heat equation, or heat transfer equation. 

For the heat equation we are interested in the change of temperature $T(x,t)$ (which is strongly related to the internal energy of a body) in space and time, hence our choice for the unknown quantity is $u = T$. Next, we need a model for the flux. Jean Babtiste Joseph Fourier (1822) gave a model for the heat flux, where the flux of heat is negative proportional to the temperature difference (heat goes from high temperatures to lower temperatures), i.e. 

\begin{align}
\vec{f}(T,x,t) &\sim 
-\Nabla T = 
\begin{bmatrix}
\partial_{x_1} T \\
\partial_{x_2} T \\
\vdots \\
\partial_{x_n} T
\end{bmatrix}
\\
\vec{f}(T) &= -d \Nabla T
\end{align}
where $d$ is the heat conduction coefficient with $d=d(T,x,t) \ge 0$ in general.

@@colbox-blue
**Remark 5:** A very simple version of the heat equation results in 1D space ($x_1 = x$) with a constant diffusion coefficient $d=const$
$$
T_t - \frac{\partial}{\partial x}(d\,T_x) = T_t - d\,T_{xx} = 0.
$$
@@

For scalar PDEs with two independent variables $x$ and $t$ of second order
(the maximum derivatives are second-order derivatives) it is common to distinguish
between different types of PDEs. Assuming a scalar second order PDE of the general form 
$$
a u_{xx} + b u_{xt} + c u_{tt} + d u_x + e u_t + f u + g = 0,
$$
we can define the quantity 
$$
\Delta(x,t) = a(x,t)\,c(x,t) - \frac{b(x,t)^2}{4}
$$
to get a classification of the different types

(i) $\Delta(x,t) > 0$, the PDE is elliptic in $(x,t)$

(ii) $\Delta(x,t) = 0$, the PDE is parabolic in $(x,t)$

(iii) $\Delta(x,t) < 0$, the PDE is hyperbolic in $(x,t)$

These definitions are motivated by the definition of conic sections. It is interesting to note, that the type of the PDE has strong implications about the solution behaviour. Hence, there are typical processes that are in nature either elliptic (e.g. gravity), hyperbolic (wave equations, advection), or parabolic (heat equation):

If we consider again our simple 1D constant coefficient heat equation and compute its type, we get $\Delta(x,t) = 0,\;\forall (x,t)$. Hence, the **heat equation is a parabolic PDE**. 

In general, parabolic PDEs model processes that evolve in time and are **not** reversible in time, such as friction, diffusion, dissipation, etc. The effect over time is a "smearing" of the quantity, smoothing out of large gradients and extrema. 

