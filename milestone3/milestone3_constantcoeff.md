+++
title = "Milestone 3 - Constant-coefficient EBM"
hascode = true
rss = "Description"
rss_title = "Milestone 3"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 3 - Constant-coefficient EBM

\toc

---

## The zero-D energy balance model

We recall that the EBM we have derived
so far is given by
$$
C(x) \partialderiv{T}{t} + A(CO_2) + B T = S_{sol}(x,t).
$$

This equation is somewhat difficult to solve analytically, because of the complexity in the solar forcing term.
To get a feeling for the analytical behaviour, as a first step, we introduce a simplification and consider a constant coefficient approximation by getting rid of the spatial and temporal dependence of the heat capacity and solar forcing coefficients: $C(x)$ and $S_{sol}(x,t) = (1 - \alpha(x)) S(x,t)$.
We do this by computing area averages in
space and averages in time to obtain
$$\label{eq:constant_coeff_ebm}
\overline{C} \partialderiv{T_A}{t} + A(CO_2) + B\, T_A = \widehat{\overline{S_{sol}}},
$$
where $\overline{C}$ is the spatial average of the heat capacity coefficient, $\widehat{\overline{S_{sol}}}$ is a spatial and temporal average of the solar forcing term, and $T_A$ is an approximation to the average of Earth's temperature.

Note that spatial averaging needs to account for the spherical shape of Earth, i.e., the spherical coordinate system. The computation of the spatial and temporal averages is detailed in [Milestone 3 - Averages](/milestone3/milestone3_averages/).

@@colbox-blue
**Remark:** Often in numerical analysis of (partial) differential equations, the important dimensions of the equations are the spatial dimensions (as these make the numerical approximation computationally expensive). Thus, one only counts the spatial dimensions to characterize the problem. We follow in the same vein and note that the simplified model only depends on time and not on any spatial coordinate after averaging. Hence, we also coin this type of model a $0D$-EBM.
@@



## Analytical solution of the $0D$-EBM

We can define the steady-state solution (also known as constant equilibrium solution) $T_{eq}$ by assuming $\partial T_A / \partial t =0$, to obtain
$$
T_{eq} = \frac{\widehat{\overline{S_{sol}}}- A}{B}.
$$

@@colbox-blue
**Remark:** This is the same temperature formula that we obtained in our simple radiation energy balance model.
@@

The ordinary differential equation \eqref{eq:constant_coeff_ebm} can be recast into
$$
\overline{C} \partialderiv{T_A}{t} = B {T_{eq}} - T_A(t),
$$
and solved analytically as
$$
T_A(t) = T_{eq} + (T_A(t=0) - T_{eq}) e^{-t/\tau},
$$
where $\tau = \overline{C}/B$ and $T_A(t=0)$ is the initial temperature of the system.

For large times, $t \rightarrow \infty$, the term $e^{-t/\tau} \rightarrow 0$, which shows that the solution will get to an equilibrium for large times.
Depending on the choice of $T_A(t=0)$, we will converge to $T_{eq}$ from below or from above:

\fig{/assets/milestone3/exp_T.png}