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

## Diffusion Term

### Modeling Physical Processes with Diffusion

* Atmospheric and oceanic currents of Earth follow some patterns, but are turbulent in nature.
* Mixing effect in turbulence.

### Diffusion Operator in Spherical Coordinates

The gradient and divergence operators are defined in Cartesian coordiates,
\begin{align}\label{eq:opsCartesian}
\Nabla T &:= \partialderiv{T}{x} \hat{x} + \partialderiv{T}{y} \hat{y} +  \partialderiv{T}{z} \hat{z},
\\
\Nabla \cdot \vec{F} &:= \partialderiv{F_x}{x} + \partialderiv{F_y}{y} + \partialderiv{F_z}{z},
\end{align}
where $T$ is a scalar field, $\hat{x}$, $\hat{y}$, and $\hat{z}$ are the unit vectors in the $x$, $y$, and $z$ directions, respectively, and $\vec{F}=F_x \hat{x} + F_y \hat{y} + F_z \hat{z} = (F_x, F_y, F_z)$ is a vector field.

We can transform the gradient and divergence operators to the spherical coordinate system by computing the relation between the unit vectors and partial derivatives of both coordinate systems.

For instance, for the unit vectors we have
<!-- \begin{align*}
\hat{r} &= \sin \colat \cos \long \hat{x} + \sin \colat \sin \long \hat{y} + \cos \colat \hat{z}
\\
\hat{\colat} &= \cos \colat \cos \long \hat x + \cos \colat \sin \long \hat{y} - \sin \colat \hat z
\\
\hat{\long} &= -\sin \long \hat x+ \cos \long \hat y
\end{align*} -->
\begin{align}\label{eq:unitvectors}
\hat{x} &= \sin \colat \cos \long \hat{r} + \cos \colat \cos \long \hat{\colat} - \sin \colat \hat{\long}
\\
\hat{y} &= \sin \colat \sin \long \hat{r} + \cos \colat \sin \long \hat{\colat} + \cos \long \hat{\long}
\\
\hat{z} &= \cos \colat \hat{r}+ \sin \colat \hat{\colat},
\end{align}
and for the partial derivatives we have
\begin{align}\label{eq:partialderivs}
\partialderiv{}{x} &= \sin \colat \cos \long \partialderiv{}{r} + \frac{\cos \colat \cos \long}{r} \partialderiv{}{\colat} + \frac{\sin \long}{r \sin \colat} \partialderiv{}{\long}
\\
\partialderiv{}{y} &= \sin \colat \sin \long \partialderiv{}{r} + \frac{\cos \colat \sin \long}{r} \partialderiv{}{\colat} + \frac{\cos \long}{r \sin \colat} \partialderiv{}{\long}
\\
\partialderiv{}{z} &= \cos \colat \partialderiv{}{r} -  \frac{\sin \colat}{r} \partialderiv{}{\colat}.
\end{align}

Inserting \eqref{eq:unitvectors} and \eqref{eq:partialderivs} into \eqref{eq:opsCartesian}, we obtain the gradient and divergence operators in spherical coordinates,
\begin{align}
\Nabla T &= \frac{\partial T}{\partial r} \hat{r} + \frac{1}{r\sin(\colat)}\frac{\partial T}{\partial \long} \hat{\long} + \frac{1}{r}\frac{\partial T}{\partial \colat} \hat{\colat}
\\
\Nabla \cdot \vec{F} &= 
\frac{1}{r^{2}}\frac{1}{\partial r}({F}_r r^{2})
 + 
\frac{1}{r\sin(\colat)}\frac{\partial}{\partial \colat}  ({F}_{\colat}\sin(\colat))
+
\frac{1}{r\sin(\colat)} \frac{\partial {F}_{\long}}{\partial \long}
\end{align}


The gradient operator is also commonly written in the vector form $(r,\long,\colat)$ as
\begin{align}
\label{eq:12}
\Nabla T = \left(\frac{\partial T}{\partial r}, \frac{1}{r\sin(\colat)}\frac{\partial T}{\partial \long}, \frac{1}{r}\frac{\partial T}{\partial \colat}\right),
\end{align}


We obtain the diffusion operator in spherical coordinates by combining the expressions for the gradient and divergence operators to obtain
\begin{align}
\Nabla \cdot (D\Nabla T) = 
\frac{1}{r^2} \left( D r^2 \partialderiv{T}{r} \right)
+
\frac{\csc^{2}(\colat)}{r^2} \frac{\partial}{\partial \long}\Bigl(D\frac{\partial T}{\partial \long}\Bigr) 
+ 
\frac{\csc (\colat)}{r^2}
\frac{\partial}{\partial \colat}\left(D \sin (\colat)\frac{\partial T}{\partial \colat}\right)
\end{align}

In our $2D$ EBM, we neglect the changes in the radial direction, define the _scaled diffusion coefficient_ $\tilde D := D / R_E^2$, and simplify the expression to obtain
\begin{align}\label{eq:diffterm}
\Nabla \cdot (D\Nabla T) = 
\underbrace{
    \csc^{2}(\colat) \frac{\partial}{\partial \long}\Bigl(\tilde D\frac{\partial T}{\partial \long}\Bigr) 
}_{\text{Term}~1}
+ 
\underbrace{
    \frac{\partial}{\partial \colat}\Bigl(\tilde D\frac{\partial T}{\partial \colat}\Bigr)
}_{\text{Term}~2}
+ 
\underbrace{
    \cot(\colat)\tilde D\frac{\partial T}{\partial \colat}.
}_{\text{Term}~3}
\end{align}

@@colbox-blue
**Attention:** The diffusion coefficient that we use in our $2D$ EBM model is scaled with Earth radius: $\tilde D := D / R_E^2$.
@@

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

## Stability Theory for Systems

* TODO: Explain mapping from $\nlat \times \nlong$ operator to a $1D$ vector. Julia command `vec` and numpy (Python) commands `.flatten()` and `reshape`.

### Computing Jacobian Matrices

#### Numerical Jacobian

It is possible to expand a general nonlinear operator, $\mathbf{R}(\mathbf{T})$, with a Taylor series to obtain,
$$ \label{eq:JacTaylorExpansion}
\mathbf{R}(\mathbf{T}) =\mathbf{R}(\mathbf{T}_0) + \mat{A(\mathbf{T}_0)} \Delta \mathbf{T}
+ \mathcal{O} ((\Delta \mathbf{T})^2),
$$
where the Jacobian matrix is defined as
$$
\mat{A} (\mathbf{T}_0) \coloneqq 
\left. \partialderiv{\mathbf{R} (\mathbf{T})}{\mathbf{T}} \right|_{\mathbf{T}_0}.
$$


This Taylor expansion can also be used to approximate the Jacobian matrix. For instance, one can assume that $\mathbf{T}$ differs from $\mathbf{T}_0$ only in one degree of freedom by a small quantity $\epsilon$, 
\begin{equation}
\mathbf{T} = \mathbf{T}_0 + \epsilon \hat{\mathbf{K}}^j,
\end{equation}
where $\hat{\mathbf{K}}$ is a vector whose entries are all zero, except in the $j$ position, where it is equal to one:
\begin{equation}
\hat{k}^j_i = 
\begin{cases} 
1 & \mathrm{if} \  i = j \\
0 & \mathrm{otherwise}.
\end{cases} 
\end{equation}

Equation \eqref{eq:JacTaylorExpansion} can then be reorganized as
\begin{equation} \label{eq:NumJacOrder1}
\mat{A} \hat{\mathbf{K}}^j = \frac{\mathbf{R}(\mathbf{T}_0+\epsilon \hat{\mathbf{K}}^j)-\mathbf{R}(\mathbf{T}_0)}{\epsilon} + \mathcal{O}(\epsilon ^ 2).
\end{equation}

For general nonlinear operators, the quantity $\mat{A} \hat{\mathbf{K}}^j$ is a good approximation for the $j^{\mathrm{th}}$ column of the Jacobian matrix if $\epsilon$ is small enough.
Since the right-hand-side term of our $2D$ EBM is a linear operator, \eqref{eq:NumJacOrder1} is exact for any $\epsilon$ and any linearization temperature $\mathbf{T}_0$, and the second and higher-order terms vanish. As a result, we can write a column of the Jacobian of our $2D$ EBM as
$$ \label{eq:NumJacEBM}
\mat{A}_{\text{EBM}} \hat{\mathbf{K}}^j = \mathbf{R}(\hat{\mathbf{K}}^j)-\mathbf{R}(\mathbf{0}).
$$

As a result, the whole Jacobian matrix can be recovered by computing \eqref{eq:NumJacEBM} $\nlat \times \nlong + 1$ times:
$$
\mat{A}_{\text{EBM}} = \left[ \mat{A}_{\text{EBM}} \hat{\mathbf{K}}^{1 \rightarrow (1,1)}, \mat{A}_{\text{EBM}} \hat{\mathbf{K}}^{2 \rightarrow (1,2)}, \ldots, \mat{A}_{\text{EBM}} \hat{\mathbf{K}}^{\nlat \times \nlong \rightarrow (\nlat, \nlong)}  \right]
$$

#### Bonus: Analytical Jacobian

