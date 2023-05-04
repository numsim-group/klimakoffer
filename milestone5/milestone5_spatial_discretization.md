+++
title = "Milestone 5"
hascode = true
rss = "Description"
rss_title = "Milestone 5"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 5

\toc

## Spatial Discretization Scheme

The challenge now is to evaluate the diffusion operator in spherical coordinates \eqref{eq:diffterm}, such that we can use it in our $2D$ EBM model.
To do that, we will obtain a _discrete_ version of the spatial derivatives in \eqref{eq:diffterm} that we can compute using the nodal values for the temperature and diffusion coefficients in our $2D$ latitude/longitude grid.

### Finite Difference Discretization of First and Second-Order Derivatives

The finite difference method is one of the simplest numerical methods to discretize derivatives. 
There are many ways to derive finite difference formulas, but we will focus on the method of the Taylor expansion, as it will give us information about the approximation errors.

@@colbox-blue
**Taylor series:** The Taylor series (or Taylor expansion) of a function is an infinite sum of polynomial terms that are expressed in terms of the function's derivatives at a point.
For instance, the Taylor deries of the function $f(x)$ around the point $x_0$ is defined as
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

**TODO: Add figure of grid with $T_{i-1}$, $T_{i}$, and $T_{i+1}$**

Using \eqref{eq:taylor}, we can write Taylor expansions for $T_{i+1}$ around $T_{i}$:
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
\deriv{T_i}{x} \approx \frac{T_{i+1} - T_{i}}{\Delta x} 
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

### Application to the Diffusion Operator in Spherical Coordinates
 
Our diffusion operator in spherical coordinates \eqref{eq:diffterm} contains a term (Term 3) of the form
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

We now have enough tools to obtain a discrete form for the diffusion term of the $2D$ EBM \eqref{eq:diffterm}.
Our goal is to obtain a finite difference scheme to discretize \eqref{eq:diffterm} with the following properties:
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

Garhering everything, and introducing $h = \Delta \lat = \Delta \long$ as the uniform mesh spacing, the central second-order finite-difference discretization of the different terms of the diffusion operator in spherical coordinates \eqref{eq:diffterm} reads
* Term 1:
\begin{align}\label{eq:disc_term1}
\left[\csc^{2}(\colat)\frac{\partial}{\partial \long}\biggl(\tilde D\frac{\partial T}{\partial \long}\Bigr)\right]_{j,i}
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
\frac{\partial}{\partial \colat}\Bigl(\tilde D\frac{\partial T}{\partial \colat}\Bigr)
\right]_{j,i}
 \approx \frac{1}{h^{2}}\biggl(-2\tilde D_{j,i}T_{j,i} 
& + \Bigl[\tilde D_{j,i} - \frac{1}{4}(\tilde D_{j+1,i} - \tilde D_{j-1,i})\Bigr] T_{j-1,i}\\
& + \Bigl[\tilde D_{j,i} + \frac{1}{4}(\tilde D_{j+1,i} - \tilde D_{j-1,i})\Bigr] T_{j+1,i}\biggr)
\end{align}

* Term 3:
\begin{align}\label{eq:disc_term3}
\left[
\cot(\colat)\tilde D\frac{\partial T}{\partial \colat}
\right]_{j,i}
 \approx \cot(\colat_{j,i})\frac{\tilde D_{j,i}}{2h}[T_{j+1,i} - T_{j-1,i}].
\end{align}

@@colbox-blue
**Remark:** We cannot use equations \eqref{eq:disc_term1}, \eqref{eq:disc_term2}, and \eqref{eq:disc_term3} for the nodes at the pole becuse of the singularity of our grid.
In fact, \eqref{eq:disc_term1} and \eqref{eq:disc_term3} are not well defined for $\colat = 0$, and the evaluation of \eqref{eq:disc_term2} is complicated because it requires a nodal value across the pole.
@@

### A Solution to the Pole Problem

To deal with the pole problem, we use the techniques proposed in the following papers:

> Bates, J. R., Semazzi, F. H. M., Higgins, R. W., & Barros, S. R. (1990). Integration of the shallow water equations on the sphere using a vector semi-Lagrangian scheme with a multigrid solver. Monthly Weather Review, 118(8), 1615-1627.

> [Barros, S. R., & Garcia, C. I. (2004). A global semi-implicit semi-Lagrangian shallow-water model on locally refined grids. Monthly weather review, 132(1), 53-65.](https://journals.ametsoc.org/downloadpdf/journals/mwre/132/1/1520-0493_2004_132_0053_agsssm_2.0.co_2.xml)

\fig{/assets/milestone5/Pole.png}

### Semi-Discrete EBM

**TODO: Write the discrete equation to solve**
