+++
title = "Milestone 3"
hascode = true
rss = "Description"
rss_title = "Milestone 3"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 3

\toc

---

@@colbox-blue
**Comment:**  In the last milestone, we focused
on the modeling aspect and introduced several
parametrizations of ,e.g., heat capacity, albedo
radiation, etc.

Before moving on to model the diffusion
operator and the diffusion coefficient we
take a step back and move deeper into
understanding the solution behaviour of the simplified
EBM we discussed in last chapter. This will also helps us
to introduce somenumerical concepts and
algorithms that we will use later forthe full
EBM.
@@

## Analytical solution of the constant-coefficient time-dependent EBM

We recall that the EBM we have derived
sofar is given by
$$
C(x) \partialderiv{T}{t} + A(CO_2) + B T = (1 - \alpha(x)) S(X,t).
$$

This equation is difficult to solve analytically.
Hence, we introduce a simplification and consider a constant coefficient approximation
by getting rid of the spacial and temporal dependence
of the heat capacity and solar forcing coefficients: $C(x)$ and $(1 - \alpha(x)) S(X,t)$.
We do this by computing the area average in
space and also averages in time.
We consider directly $C(x)$ and $(1 - \alpha(x)) S(X,t)$
as computed in milestone 2.

## Area average

We have to be careful about computing the
area average: while we use a Cartesian computational grid in longitude and latitude, we need the average over the
sphere surface.

We consider the latitude/longitude grid:

\fig{/assets/milestone3/Grid.png}

In the figure, the discretization nodes are marked in purple, the grid _lines_ are marked in gray, and we have drawn additional auxiliary grid lines (dotted lines in red) between our computational nodes.
We are going to use these auxiliary grid lines to define the surface area that corresponds to each computational node.

Observe that, even though we use a regular grid, the area that corresponds to each computational node of our grid depends only on the latitude, and not on the longitude.
Therefore, our goal is to compute the corresponding
approximation of the area for each latitude
grid coordinate ($\lat_j$).

In addition, observe that each node at the pole maps to a triangle, and not to a quadrilateral.
In preparation for the discretization of the diffusion operator in sperical
coordinates ([milestone 5](/milestone5/milestone5_menu/)), we already start to
treat the pole regions ($j=1$ and $j=\nlat$)
in a special way.

### Poles (North and South)

We need to compute the area of a spherical cap at the angle $\lat = \frac{\Delta \lat}{2}$.

\fig{/assets/milestone3/SphericalCap.png}
* Source: Wikipedia

Using mathematical analysis, we can derive the area of a spherical cap with radius $R_E$ as
$$
\tilde A_{\text{pole}} = 2 \pi R_E^2 (1-\cos \lat) = 2 \pi R_E^2 \left(1-\cos \left(\frac{\Delta \lat}{2}\right)\right).
$$

As we will use the area calculation only in to compute area averages of quantities on the surface of the spghere, we will directly normalize this value with the surface area of the Earth (spehere), i.e.
$$\label{eq:cap}
A_{\text{pole}} = \frac{1}{4\pi R_E^2} \tilde A_{\text{pole}} = \frac{1}{2} \left(1-\cos \left(\frac{\Delta \lat}{2}\right)\right).
$$

### Interior nodes

For the interior nodes, we consider a spherical segment between the two angles $\lat_j - \Delta \lat/2$ and $\lat_j + \Delta \lat/2$.
We can use \eqref{eq:cap} to compute the difference between two spherical caps and obtain
\begin{align}
A_{\text{int}} &= \frac{1}{2} \left(\cos \left(\lat_j - \frac{\Delta \lat}{2}\right) - \cos \left(\lat_j + \frac{\Delta \lat}{2}\right)\right)
\\
&= \frac{1}{2}  \sin \left(\frac{\Delta \lat}{2}\right) \sin \left(\lat_j\right).
\end{align}