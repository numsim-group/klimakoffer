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

## Planetary Heat Transfer

### Physical Processes in the Atmosphere an Oceans

* Earth's climate is basically the transfer of heat from the (hot) equator to the (cold) poles.
* Atmospheric currents of Earth follow some patterns: Hadley, Ferrel, and polar cells. 
* Effect of Coreolis force in large-scale vortices(?)
* Atmospheric currents are turbulent in nature.
* Oceanic currents can be classified in surface and deep-ocean currents.

\fig{/assets/milestone5/heat_transfer_north.png}
* Heat transfer in the northern hemisphere. From "Maslin (2013). 'Climate: A Very Short Introduction'"

### Modeling Physical Processes with Diffusion

* Fluid dynamics is very complex.
* Turbulence is a multi-scale phenomenon that requires resolution of large and small eddies.
* Unfortunately, a full general recirculation model is out of the scope of the course.
* Modeling: Mixing effect in turbulence.

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
