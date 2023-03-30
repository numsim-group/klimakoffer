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
C(x) \partialderiv{T}{t} + A(CO_2) + B T = (1 - \alpha(x)) S(x,t).
$$

This equation is difficult to solve analytically.
Hence, we introduce a simplification and consider a constant coefficient approximation
by getting rid of the spatial and temporal dependence
of the heat capacity and solar forcing coefficients: $C(x)$ and $(1 - \alpha(x)) S(X,t)$.
We do this by computing area averages in
space and averages in time to obtain
$$\label{eq:constant_coeff_ebm}
\overline{C} \partialderiv{\overline{T}}{t} + A(CO_2) + B \overline{T} = \widehat{\overline{Q_{\alpha} S}},
$$
where $\overline{C}$ is the spatial average of the heat capacity coefficient, $\widehat{\overline{Q_{\alpha} S}}$ is a spatial and temporal average of the solar forcing term, and $\overline{T}$ is an approximation to the spatial average of Earth's temperature.

The computation of the spatial and temporal averages is detailed in [Milestone 3 - Averages](/milestone3/milestone3_averages/).

As this model only depends on time and not on any spatial coordinate, we also coin this type of model $0D$-EBM.

## Analytical solution

We can define the steady-state solution (also known as equilibrium solution) $T_{eq}$ by assuming $\partial \overline{T} / \partial t =0$, to obtain
$$
T_{eq} = \frac{\widehat{\overline{Q_{\alpha} S}} - A}{B}.
$$

@@colbox-blue
**Remark:** This is the same temperature that we obtained in our simple radiation energy balance model.
@@

The ordinary differential equation \eqref{eq:constant_coeff_ebm} can be recast into
$$
\overline{C} \partialderiv{\overline{T}}{t} = B {T_{eq}} - \overline{T}(t),
$$
and solved analytically as
$$
\overline{T}(t) = T_{eq} + (\overline{T}_0 - T_{eq}) e^{-t/\tau},
$$
where $\tau = \overline{C}/B$ and $\overline{T}_0 \coloneqq \overline{T}(t=0)$ is the initial temperature of the system.

For large times, $t \rightarrow \infty$, the term $e^{-t/\tau} \rightarrow 0$, which shows that the solution will get to an equilibrium for large times.
Depending on the choice of $\overline{T}_0$, we will converge to $T_{eq}$ from below or from above:

\fig{/assets/milestone3/exp_T.png}