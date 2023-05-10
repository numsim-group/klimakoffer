+++
title = "Milestone 5"
hascode = true
rss = "Description"
rss_title = "Milestone 5"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 5 - Spatial Discretization Scheme

\toc

The challenge now is to evaluate the [diffusion operator in spherical coordinates](/milestone5/milestone5_heat_transfer/#eqdiffterm), such that we can use it in our $2D$ EBM model.
To do that, we will obtain a _discrete_ version of the spatial derivatives that we can compute using the nodal values for the temperature and diffusion coefficients in our $2D$ latitude/longitude grid.

## Finite Difference Discretization of First and Second-Order Derivatives

The finite difference method is one of the simplest numerical methods to discretize derivatives. 
There are many ways to derive finite difference formulas, but we will focus on the method of the Taylor expansion, as it will give us information about the approximation errors.

@@colbox-blue
**Taylor series:** The Taylor series (or Taylor expansion) of a function is an infinite sum of polynomial terms that are expressed in terms of the function's derivatives at a point.
For instance, the Taylor series of the function $f(x)$ around the point $x_0$ is defined as
\begin{align}\label{eq:taylor}
f(x) 
&= f(x_0) + \deriv{f(x_0)}{x} (x-x_0) + \frac{1}{2!} \deriv{^2f(x_0)}{x^2} (x-x_0)^2 + \frac{1}{3!}\deriv{^3f(x_0)}{x^3} (x-x_0)^3 + \cdots
\\ 
&= f(x_0) + \deriv{f(x_0)}{x} \Delta x + \frac{1}{2!} \deriv{^2f(x_0)}{x^2} \Delta x^2 + \frac{1}{2!} \deriv{^3f(x_0)}{x^3} \Delta x^3 + \cdots
\\
&= \sum_{n=1}^{\infty} \frac{f^{(n)}(x_0)}{n!} \Delta x,
\end{align}
where $f^{(n)}(x_0)$ is a short-hand notation for the $n^{\text{th}}$ derivative of function $f(x)$ at $x_0$, and the distance between the evaluation point and the point $x$ and $x_0$ is $\Delta x \coloneqq x - x_0$.
@@

Let us consider a uniform one-dimensional grid containing temperature values.

\fig{/assets/milestone5/spatial_disc_1d.png}

Using \eqref{eq:taylor}, we can write Taylor expansion for $T_{i+1}$ around $T_{i}$:
$$\label{eq:taylortemp1}
T_{i+1} = T_i + \deriv{T_i}{x} \Delta x + \frac{1}{2} \deriv{^2T_i}{x^2} \Delta x^2 + \frac{1}{6} \deriv{^3T_i}{x^3} \Delta x^3 + \cdots.
$$

By manipulating \eqref{eq:taylortemp1}, we can obtain an expression for the first-order derivative,
\begin{align}
\deriv{T_i}{x} &= \frac{T_{i+1} - T_{i}}{\Delta x} - \frac{1}{2} \deriv{^2T_i}{x^2} \Delta x + \ldots
\\
&= \frac{T_{i+1} - T_{i}}{\Delta x} + \mathcal{O}(\Delta x),
\end{align}
where $\mathcal{O}(\Delta x)$ indicates that the largest term is of the order of $\Delta x$. 
Hence, if we truncated the Taylor series there, 
$$
\deriv{T_i}{x} \approx \frac{T_{i+1} - T_{i}}{\Delta x},
$$
we would obtain an approximation of the derivative $\d T_i / \d x$, which would converge to the exact value of the derivative with first-order accuracy (i.e., $\text{error} \sim \Delta x^1$) as the grid size is reduced, $\Delta x \rightarrow 0$.

It is also possible to write Taylor expansions of $T_{i-1}$ around $T_{i}$,
$$\label{eq:taylortemp2}
T_{i-1} = T_i - \deriv{T_i}{x} \Delta x + \frac{1}{2} \deriv{^2T_i}{x^2} \Delta x^2 - \frac{1}{6} \deriv{^3T_i}{x^3} \Delta x^3 + \cdots,
$$
or $T_{i+2}$ around $T_{i}$,
$$\label{eq:taylortemp3}
T_{i+2} = T_i + 2 \deriv{T_i}{x} \Delta x + 2 \deriv{^2T_i}{x^2} \Delta x^2 + \frac{4}{3} \deriv{^3T_i}{x^3} \Delta x^3 + \cdots,
$$
etc.

Using \eqref{eq:taylortemp1}, \eqref{eq:taylortemp2}, and \eqref{eq:taylortemp3}, we can obtain the following approximations for the **first derivative** at $x_i$:
* Right-sided finite-difference scheme of first order ($\text{error} \sim \Delta x^1$):
$$
\deriv{T_i}{x} \approx \frac{T_{i+1} - T_{i}}{\Delta x}.
$$
* Left-sided finite-difference scheme of first order ($\text{error} \sim \Delta x^1$):
$$
\deriv{T_i}{x} \approx \frac{T_{i} - T_{i-1}}{\Delta x}.
$$
* Central finite-difference scheme of second order ($\text{error} \sim \Delta x^2$):
$$\label{eq:central_first}
\deriv{T_i}{x} \approx \frac{T_{i+1} - T_{i-1}}{2 \Delta x}.
$$
* Right-sided finite-difference scheme of second order ($\text{error} \sim \Delta x^2$):
$$
\deriv{T_i}{x} \approx \frac{-3 T_i + 4T_{i+1} - T_{i+2}}{2 \Delta x}.
$$
* And many other possibilities.

@@colbox-blue
**Remark:** High-order approximations are desired as the higher the order, the faster we approach the exact value of the derivative as we refine the grid. However, high-order accuracy sometimes come with the price of more operations to compute.
@@

@@colbox-blue
**Remark:** Note that some finite difference approximations are central, and some are left-sided or right-sided. Depending on how the finite-difference formula is defined, it will have some bias as to how information travels in the domain.
@@

Using \eqref{eq:taylortemp1}, \eqref{eq:taylortemp2}, and \eqref{eq:taylortemp3}, we can also obtain approximations for the **second derivative** at $x_i$:
* Right-sided finite-difference scheme of first order ($\text{error} \sim \Delta x^1$):
$$
\deriv{^2T_i}{x^2} \approx \frac{T_{i+2} - 2T_{i+1} + T_{i}}{\Delta x^2}.
$$
* Central finite-difference scheme of second order ($\text{error} \sim \Delta x^2$):
$$\label{eq:central_second}
\deriv{^2T_i}{x^2} \approx \frac{T_{i+1} - 2T_i + T_{i-1}}{\Delta x^2}.
$$
* And many other possibilities.

## Application to the Diffusion Operator in Spherical Coordinates
 
Our [diffusion operator in spherical coordinates](/milestone5/milestone5_heat_transfer/#eqdiffterm) contains a term (Term 3) of the form
$$
\mathcal{L}_1(T) = \partialderiv{T}{x},
$$
which can be discretized with the first-order derivative finite-difference approximations, and two terms (Terms 1 and 2) of the form
$$
\mathcal{L}_2(T) = \partialderiv{}{x} \left( \tilde D \partialderiv{T}{x} \right),
$$
which is slightly more involved.

One possibility to discretize $\mathcal{L}_2(T)$ is to apply the chain rule of differentiation at the continuous level,
$$
\mathcal{L}_2(T) = \partialderiv{}{x} \left( \tilde D \partialderiv{T}{x} \right) = 
\partialderiv{\tilde D}{x} \partialderiv{T}{x} +
\tilde D \partialderiv{^2T}{x^2},
$$
such that it is expressed in terms of first and second-order derivatives.

@@colbox-blue
**Remark:** Even though we are now dealing with partial derivatives, the finite difference formulas can be applied dimension by dimension.
@@

We now have enough tools to obtain a discrete form for the [diffusion term of the $2D$ EBM](/milestone5/milestone5_heat_transfer/#eqdiffterm).
Our goal is to obtain a finite difference scheme with the following properties:
* The scheme must be **central**: To model the parabolic nature of the heat conduction term (see [Milestone 2 - Conservation and Balance](/milestone2/milestone2_conservation-and-balance/)), we require the scheme to be symmetric. 
* The scheme must be **compact**: To avoid wide stencils, we require the scheme to only use the information from neighboring nodes.
* The scheme must be as **accurate** as possible: We will select a scheme of the highest order possible, such that our approximation is as accurate as possible given our constraints.

Given these criteria, we select the central second-order accurate finite-difference approximations for the first-order \eqref{eq:central_first} and second-order \eqref{eq:central_second} derivatives, to obtain for the degree of freedom $i$
$$
\mathcal{L}_1(T_i) = \tilde D_i \partialderiv{T_i}{x} \approx \tilde D_i \frac{T_{i+1} - T_{i-1}}{2 \Delta x},
$$
and
\begin{align}
\mathcal{L}_2(T_i) = \partialderiv{}{x} \left( \tilde D_i \partialderiv{T_i}{x} \right) =& 
\partialderiv{\tilde D}{x} \partialderiv{T}{x} +
\tilde D \partialderiv{^2T}{x^2}
\\
\approx &
\frac{\tilde D_{i+1} - \tilde D_{i-1}}{2 \Delta x}
\frac{T_{i+1} - T_{i-1}}{2 \Delta x}
+
\tilde D_i \frac{T_{i+1} - 2T_i + T_{i-1}}{\Delta x^2},
\end{align}
which simplifies to
\begin{align}
\mathcal{L}_2(T_i) = 
\approx 
\frac{1}{\Delta x^2}
\biggl(
-2 \tilde D_i T_i 
&+ \Bigl[\tilde D_i + \frac{1}{4} \left( \tilde D_{i+1} - \tilde D_{i-1} \right)\Bigr] T_{i+1}
\\
&+ 
\Bigl[\tilde D_i - \frac{1}{4} \left( \tilde D_{i+1} - \tilde D_{i-1} \right)\Bigr] T_{i-1}
\biggr).
\end{align}

Garhering everything, and introducing $h = \Delta \lat = \Delta \long$ as the uniform mesh spacing, the central second-order finite-difference discretization of the different terms of the [diffusion operator in spherical coordinates](/milestone5/milestone5_heat_transfer/#eqdiffterm) reads
$$\label{eq:diffop_innernodes}
\diffop_{j,i} = 
\underbrace{
\left[\csc^{2}(\colat)\frac{\partial}{\partial \long}\biggl(\tilde D (\colat,\long)\frac{\partial T}{\partial \long}\Bigr)\right]_{j,i}
}_{\text{Term~1}}
+
\underbrace{
\left[
\frac{\partial}{\partial \colat}\Bigl(\tilde D (\colat,\long)\frac{\partial T}{\partial \colat}\Bigr)
\right]_{j,i}
}_{\text{Term~2}}
+
\underbrace{
\left[
\cot(\colat)\tilde D (\colat,\long)\frac{\partial T}{\partial \colat}
\right]_{j,i}
}_{\text{Term~3}},
$$
with the terms
* Term 1:
\begin{align}\label{eq:disc_term1}
\left[\csc^{2}(\colat)\frac{\partial}{\partial \long}\biggl(\tilde D (\colat,\long)\frac{\partial T}{\partial \long}\Bigr)\right]_{j,i}
\approx 
\frac{\csc^{2}(\colat_{j,i})}{h^2}
\biggl(-2\tilde D_{j,i}T_{j,i}
&
+  \Bigl[\tilde D_{j,i} - \frac{1}{4}(\tilde D_{j,i+1} - \tilde D_{j,i-1})\Bigr]T_{j,i-1} \\
& 
+ \Bigl[\tilde D_{j,i} + \frac{1}{4}(\tilde D_{j,i+1} - \tilde D_{j,i-1})\Bigr] T_{j,i+1}\biggr)
\end{align}

* Term 2:
\begin{align}\label{eq:disc_term2}
\left[
\frac{\partial}{\partial \colat}\Bigl(\tilde D (\colat,\long)\frac{\partial T}{\partial \colat}\Bigr)
\right]_{j,i}
 \approx \frac{1}{h^{2}}\biggl(-2\tilde D_{j,i}T_{j,i} 
& + \Bigl[\tilde D_{j,i} - \frac{1}{4}(\tilde D_{j+1,i} - \tilde D_{j-1,i})\Bigr] T_{j-1,i}\\
& + \Bigl[\tilde D_{j,i} + \frac{1}{4}(\tilde D_{j+1,i} - \tilde D_{j-1,i})\Bigr] T_{j+1,i}\biggr)
\end{align}

* Term 3:
\begin{align}\label{eq:disc_term3}
\left[
\cot(\colat)\tilde D (\colat,\long)\frac{\partial T}{\partial \colat}
\right]_{j,i}
 \approx \cot(\colat_{j,i})\frac{\tilde D_{j,i}}{2h}[T_{j+1,i} - T_{j-1,i}].
\end{align}

@@colbox-blue
**Remark:** We cannot use equations \eqref{eq:disc_term1}, \eqref{eq:disc_term2}, and \eqref{eq:disc_term3} for the nodes at the pole becuse of the singularity of our grid.
In fact, \eqref{eq:disc_term1} and \eqref{eq:disc_term3} are not well defined for $\colat = 0$, and the evaluation of \eqref{eq:disc_term2} is complicated because it requires a nodal value across the pole.
@@

## A Solution to the Pole Problem

To deal with the pole problem, we use the techniques proposed in the following papers:

> Bates, J. R., Semazzi, F. H. M., Higgins, R. W., & Barros, S. R. (1990). Integration of the shallow water equations on the sphere using a vector semi-Lagrangian scheme with a multigrid solver. Monthly Weather Review, 118(8), 1615-1627.

> [Barros, S. R., & Garcia, C. I. (2004). A global semi-implicit semi-Lagrangian shallow-water model on locally refined grids. Monthly weather review, 132(1), 53-65.](https://journals.ametsoc.org/downloadpdf/journals/mwre/132/1/1520-0493_2004_132_0053_agsssm_2.0.co_2.xml)

We start by considering the [diffusion operator in spherical coordinates](/milestone5/milestone5_heat_transfer/#eqdiffterm), but we combine Terms 2 and 3 to obtain
\begin{align}\label{eq:difftermpole}
\diffop(\colat,\long) \coloneqq 
\Nabla \cdot (D\Nabla T) = 
    \csc^{2}(\colat) \frac{\partial}{\partial \long}\Bigl(\tilde D (\colat,\long)\frac{\partial T}{\partial \long}\Bigr) 
    +
    \csc(\colat) \frac{\partial }{\partial \colat} \left( \tilde D (\colat,\long) \sin (\colat) \partialderiv{T}{\colat} \right).
\end{align}
<!-- where the newly defined variable $\diffop$ contains the rest of the terms of the EBM, i.e., the heat capacity term, the outgoinglongwave ratiation, and the solar forcing. -->

We will proceed to derive a discretization for the north pole. The discretization for the south pole can be derived in an analogous manner.
Let us consider a _control area_ of the size of the northern polar cell (green in the figure), which spans $\colat \in [0,h/2]$, $\long \in [0, 2\pi]$.

\fig{/assets/milestone5/Pole.png}

To avoid the singularity at the pole, we will now consider an integral form of \eqref{eq:difftermpole}, which we obtain by integrating the equation over the polar cell,
\begin{align}\label{eq:difftermpole_weak}
\int_0^{\frac{h}{2}} \int_0^{2\pi}\diffop(\colat,\long)R_E^2 \sin \colat \d \long \d \colat =&
    \int_0^{\frac{h}{2}}  \int_0^{2\pi}
    \csc(\colat) \frac{\partial}{\partial \long}\Bigl(\tilde D (\colat,\long)\frac{\partial T}{\partial \long}\Bigr) 
    R_E^2  \d \long \d \colat
    \\
    &+
    \int_0^{\frac{h}{2}}  \int_0^{2\pi}
     \frac{\partial }{\partial \colat} \left( \tilde D (\colat,\long) \sin (\colat) \partialderiv{T}{\colat} \right)
    R_E^2  \d \long \d \colat.
\end{align}

@@colbox-blue
**Remark:** In \eqref{eq:difftermpole_weak}, we have used the fact that the differential of area in spherical coordinates is defined as
$$
\d S = 
\left| \partialderiv{\mathbf{r}}{\colat} \times \partialderiv{\mathbf{r}}{\varphi}  \right| \d \long \d \colat
=
R_E^2 \sin \colat \d \long \d \colat,
$$
where $\mathbf{r}$ is a vector that describes the surface of Earth, i.e.,
\begin{align}
\mathbf{r} = (\hat{x}, \hat{y}, \hat{z})\bigg\vert_{r=R_E}.
\end{align}

To compute a surface integral, we use
\begin{align}
\int_S \diffop(\colat,\long) \d S &= \int_0^{\Delta \theta/2} \int_{-\pi}^{\pi}
\diffop(\colat,\long) \left| \partialderiv{\mathbf{r}}{\colat} \times \partialderiv{\mathbf{r}}{\varphi}  \right|  \d \varphi \d \colat
\\
&= \int_0^{\Delta \theta/2} \int_{-\pi}^{\pi}
\diffop(\colat,\long)
R_E^2 \sin \colat \d \varphi \d \colat.
\end{align}
@@

We can show that the first integral in the right-hand side of \eqref{eq:difftermpole_weak} vanishes using the fundamental theorem of calculus:
\begin{align}\label{eq:intform1}
\int_0^{\frac{h}{2}}  \int_0^{2\pi}
    \csc(\colat) \frac{\partial}{\partial \long}\Bigl(\tilde D (\colat,\long)\frac{\partial T}{\partial \long}\Bigr) 
    R_E^2  \d \long \d \colat
=&
R_E^2 \int_0^{\frac{h}{2}} \csc(\colat)
    \int_0^{2\pi}
    \frac{\partial}{\partial \long}\Bigl(\tilde D (\colat,\long)\frac{\partial T}{\partial \long}\Bigr) \d \long 
    \d \colat
\\
=&
R_E^2 \int_0^{\frac{h}{2}} \csc(\colat)
    \underbrace{\left[
    \tilde D (\colat,\long)\frac{\partial T}{\partial \long} \right]_0^{2\pi}}_{=0} 
    \d \colat
\\
=& 0.
\end{align}

@@colbox-blue
**Remark:** As mentioned above, there is a singularity at $\colat = 0$. As shown in \eqref{eq:intform1}, we can deal with the singularity using the integral form of the equation. As can be derived with standard calculus, the value of the integral is well defined.
@@

Similarly, we can use the fundamental theorem of calculus to integrate the second term in the right-hand side of \eqref{eq:difftermpole_weak} over the colatitude:
\begin{align}
\int_0^{\frac{h}{2}} \int_0^{2\pi} \diffop(\colat,\long) R_E^2 \sin \colat \d \long \d \colat =&
    R_E^2 \int_0^{2\pi} \int_0^{\frac{h}{2}} 
     \frac{\partial }{\partial \colat} \left( \tilde D (\colat,\long) \sin (\colat) \partialderiv{T}{\colat} \right)
     \d \colat \d \long
\\
=& 
R_E^2
\int_0^{2\pi} 
\left[  \tilde D (\colat,\long) \sin (\colat) \partialderiv{T}{\colat} \right]_0^{\frac{h}{2}} 
\d \long
\\
=& 
R_E^2
\int_0^{2\pi} 
\tilde D \left( \frac{h}{2} ,\long \right) \sin \left( \frac{h}{2} \right) \partialderiv{T}{\colat} \bigg\rvert_{\colat=\frac{h}{2}}  
\d \long.
\end{align}

As a next step, we make the assumption that $\diffop(\colat,\long)$ is a constant in the polar cell, such that we can approximate the integral on the left-hand side with a mid-point rule,
\begin{align}
\int_0^{\frac{h}{2}} \int_0^{2\pi} \diffop(\colat,\long) R_E^2 \sin \colat \d \long \d \colat
\approx&
\int_0^{\frac{h}{2}} \int_0^{2\pi} \diffop_{1,i} R_E^2 \sin \colat \d \long \d \colat
\\
=&
\diffop_{1,i}
\int_0^{\frac{h}{2}} \int_0^{2\pi}  R_E^2 \sin \colat \d \long \d \colat
\\
=&
\diffop_{1,i}
\left( 2 \pi R_E^2 \left(1-\cos \left(\frac{\Delta \lat}{2}\right) \right) \right)
\\
=&
\diffop_{1,i} (\texttt{area[1]}) (4\pi R_E^2)
\end{align}
where $\diffop_{1,i}$ is the nodal value of $\diffop(\colat,\long)$ at $\colat_j=0$ ($j=1$) for any longitude position $i$. This is an approximation of the integral that assumes that the nodal value of the mid-point corresponds to the area-average of $\diffop(\colat,\long)$ in the polar cell.
Recall that the first entry of the array $\texttt{area}$ contains the total area of the polar cap normalized with the area of the sphere (see [Milestone 3 - Averages](/milestone3/milestone3_averages/)).

We then obtain
\begin{align}\label{eq:difftermpole_weak2}
\diffop_{1,i} (\texttt{area[1]}) (4\pi R_E^2)
=& 
R_E^2
\int_0^{2\pi} 
\tilde D \left( \frac{h}{2},\long \right) \sin \left( \frac{h}{2} \right) \partialderiv{T}{\colat} \bigg\rvert_{\colat=\frac{h}{2}}  
\d \long.
\end{align}

Now, we approximate the derivative of temperature with respect to colatitude at each discrete longitude, $\long_i$, with a central finite-difference approximation of second-order accuracy,
$$
\left[ 
    \partialderiv{T}{\colat} 
\right]_{\colat=\frac{h}{2},i} \approx 
\frac{T_{2,i}-T_{1,i}}{h},
$$
and the diffusion coefficient at $\colat=\frac{h}{2}$ with an area-weighted average,
$$
\left[ \tilde D (\colat,\long)\right]_{\colat=\frac{h}{2},i} \approx
\bar{\tilde{D}}_i^{\text{NP}} =
\frac{\frac{\texttt{area[1]}}{\nlong} \tilde D_{1,i} + \texttt{area[2]}\tilde D_{2,i}}{\frac{\texttt{area[1]}}{\nlong} + \texttt{area[2]}}.
$$

Finally, we approximate the integral on the right-hand side of \eqref{eq:difftermpole_weak2} with a rectangular quadrature rule to obtain
$$\label{eq:diffop_NP}
\diffop_{1,i} 
=
\frac{\sin \left( h/2 \right)}{4\pi \, \texttt{area[1]}}
\sum_{i=1}^{\nlong}
\bar{\tilde{D}}_i^{\text{NP}}
\left[ T_{2,i}-T_{1,i} \right].
$$

Following a similar procedure, we obtain for the south pole
$$\label{eq:diffop_SP}
\diffop_{\nlat,i} 
=
\frac{\sin \left( h/2 \right)}{4\pi \, \texttt{area[\nlat]}}
\sum_{i=1}^{\nlong}
\bar{\tilde{D}}_i^{\text{SP}}
\left[ T_{\nlat-1,i}-T_{\nlat,i} \right],
$$
with
$$
\bar{\tilde{D}}_i^{\text{SP}} =
\frac{\frac{\texttt{area[\nlat]}}{\nlong} \tilde D_{\nlat,i} + \texttt{area[\nlat-1]}\tilde D_{\nlat-1,i}}{\frac{\texttt{area[\nlat]}}{\nlong} + \texttt{area[\nlat-1]}}.
$$

## Semi-Discrete EBM

With the spatial discretization that we derived above, we can rewrite the PDE that describes our [$2D$ EBM in spherical coordinates](/milestone5/milestone5_heat_transfer/#eqebm_spherical) as an ODE for each degree of freedom $(j,i)$:
$$\label{eq:semidisc_ebm_point}
\deriv{T_{j,i}}{t} = 
\underbrace{
    \frac{L_{j,i}}{C_{j,i}} - \frac{B}{C_{j,i}} T_{j,i}
 }_{R_{j,i}(T)}
+
\underbrace{
    \frac{1}{C_{j,i}} \left(S_{sol, j, i}(t) - A\right)
}_{F_{j,i}},
$$
where we have gathered the temperature-dependent terms under $R_{j,i}(T)$ and the sources and sinks under $F_{j,i}$.

In equation \eqref{eq:semidisc_ebm_point}, $C_{j,i}$ is the point-wise heat capacity that depends on the geography (see [Milestone 2 - Heat Capacity](/milestone2/milestone2_heat-capacity/)); $L_{j,i}$ is the discrete point-wise heat transfer term computed with \eqref{eq:diffop_innernodes} for the interior nodes, with \eqref{eq:diffop_NP} for the north pole ($j=1$), and with \eqref{eq:diffop_SP} for the south pole ($j=\nlat$); $A$ and $B$ are the outgoing longwave radiation parameters (see [Milestone 2 - Radiation Modeling](/milestone2/milestone2_radiation/)); and $S_{sol, j, i} = (1 - \alpha_{j,i}) S_{j,i}$ is the point-wise time-dependent solar forcing term that depends on the [point-wise insolation $S_{j,i}$](/milestone2/milestone2_solar-forcing/#eqinsolation) and the [point-wise albedo $\alpha_{j,i}$](/milestone2/milestone2_albedo/).

We can rewrite \eqref{eq:semidisc_ebm_point} in matrix form as
$$\label{eq:semidisc_ebm_matrix}
\dot{\mat{T}} = 
\mat{R}(\mat{T})
+
\mat{F}(t),
$$
where $\mat{T} \in \R^{\nlat \times \nlong}$ is a matrix that contains all point-wise values of temperature, $\mat{R}(\mat{T}) \in \R^{\nlat \times \nlong}$ is a matrix operator applied on the temperature matrix, $\mat{F}(t) \in \R^{\nlat \times \nlong}$ is the time-dependent matrix containing the sources and sinks of the method, and we use the dot operator to denote the time derivative term.