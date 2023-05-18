+++
title = "Milestone 5"
hascode = true
rss = "Description"
rss_title = "Milestone 5"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++


# Milestone 5 - Constructing a System Matrix with FD Approach

## Constructing a System Matrix with FD Approach

[(Download description in PDF format)](/assets/milestone5/description.pdf)\\
[(Download *The_World128x65.dat*)](/assets/milestone5/input/The_World128x65.dat)\\

In this milestone, you will build the heart of the climate model: the system matrix for solving the energy balance model.
The EBM model reads
\begin{equation}
C(x)\frac{\partial T}{\partial t}+A(CO_2)+BT - \underbrace{\vec{\nabla} \cdot (D(x)\nabla T)}_{L} = S_\text{sol}(x,t)
\label{EBM}
\end{equation}
in spherical coordinates on the sphere's surface with the pair $ (\tilde{\theta},\varphi)$
\begin{equation}
L(\tilde{\theta},\varphi) = \csc^2(\tilde{\theta})\frac{\partial}{\partial \varphi}\left(\tilde{D} \frac{\partial T}{\partial \varphi}\right)+\frac{\partial}{\partial \tilde{\theta}}\left(\tilde{D} \frac{\partial T}{\partial \tilde{\theta}}\right)+ \cot(\tilde{\theta})\tilde{D}\frac{\partial T}{\partial \tilde{\theta}}
\label{Operator}
\end{equation}
When we apply the finite difference method, we have for each degree of freedom:
\begin{equation}
\frac{\partial T_{j,i}}{\partial t} \approx \underbrace{\frac{L_{j,i}}{C_{j,i}} - \frac{B}{C_{j,i}} T_{j,i}}_{R_{j,i}(T)} + \underbrace{ \frac{1}{C_{j,i}} \left(S_{\text{sol}, j, i}(t) - A\right)}_{F_{j,i}}
\label{DoF}
\end{equation}
The discretization of the heat conduction term on the surface of a sphere with a Cartesian mesh of cell size $h$ (uniform for latitude and longitude) reads:

* For the inner degrees of freedom:
\begin{equation} \label{innerDoF}
\begin{split}
L_{j,i} = \frac{\csc^2(\tilde{\theta}_{j,i})}{h^2}\Big(-2\tilde{D}_{j,i}T_{j,i}&+\left[\tilde{D}_{j,i}-\frac{1}{4}(\tilde{D}_{j,i+1}-\tilde{D}_{j,i-1})\right]T_{j,i-1}  \\
& \left. + \left[\tilde{D}_{j,i}+\frac{1}{4}(\tilde{D}_{j,i+1}-\tilde{D}_{j,i-1})\right]T_{j,i+1}\right) \\
+ \frac{1}{h^2}\Big(-2\tilde{D}_{j,i}T_{j,i}&+\left[\tilde{D}_{j,i}-\frac{1}{4}(\tilde{D}_{j+1,i}-\tilde{D}_{j-1,i})\right]T_{j-1,i}  \\
& \left. +\left[\tilde{D}_{j,i}+\frac{1}{4}(\tilde{D}_{j+1,i}-\tilde{D}_{j-1,i})\right]T_{j+1,i}\right) \\
 + \cot(\tilde{\theta}_{j,i})\frac{\tilde{D}_{j,i}}{2h}&\left[T_{j+1,i}-T_{j-1,i}\right]
\end{split}
\end{equation}

*  For the degrees of freedom at the poles:
\begin{align*}
L_{1,i} = \frac{\sin(\frac{h}{2})}{4 \pi \text{area[1]}}\sum^{\text{n\_longitude}}_{k=1}\frac{1}{2}\left(\tilde{D}_{1,k}+\tilde{D}_{2,k}\right)\left[T_{2,k}-T_{1,k}\right] \\
L_{\text{n\_lat},i} = \frac{\sin(\frac{h}{2})}{4 \pi \text{area[n\_lat]}}\sum^{\text{n\_longitude}}_{k=1}\frac{1}{2}\left(\tilde{D}_{\text{n\_lat},k}+\tilde{D}_{\text{n\_lat-1},k}\right) \left[T_{\text{n\_lat-1},k}-T_{\text{n\_lat},k}\right]
\label{polesDoF}
\end{align*}
for all $i \in \{1,\dots,\text{n\_longitude}\}$ where $\texttt{area[1]} := \texttt{area[n\_lat]} := 0.5\left(1-\cos\left(\frac{h}{2}\right)\right)$ contains the normalized area of the polar cap.
For space reasons the $n\_latitude$ from the script is shortened to $n\_lat$.


To implement the model, proceed as follows:

1. Either implement the mesh as a struct in Julia or a class in Python with a constructor or alteratively write a function returning all relevant mesh data given the geography data as input.

| **Parameter** | **Description** |
|-------------|-----------|
| n\_latitude | Number of DOFs in longitude |
| n\_longitude | Number of DOFs in latitude |
| NDOF | Total number of DOFs |
| h | Grid size in radians |
| geom | Geometrical parameter at poles with $geom = \sin(\frac{h}{2})/(4\pi\text{ }area[1])$ |
| csc2 | Metric term for the transformation to spherical coordinates, the vector with entries $\csc^2(\theta_j)$|
| cot | Metric term for the transformation to spherical coordinates, the vector with entries $\cot(\theta_j)$ |
| area | Area of each cell on the surface of the sphere which only depends on latitude |

2. Calculate the earth surface type dependent heat conduction coefficients $\tilde{D}(\tilde{\theta},\varphi) \in \mathbb{R}^{n_y \times n_x}$ by implementing a function *calc\_diffusion\_coefficients* that maps the following equation:

\begin{align}
\tilde{D}(\tilde{\theta},\varphi) \in \mathbb{R}^{n_y \times n_x} =
\begin{cases}
\tilde{D}_{ocean,poles}+ (\tilde{D}_{ocean,equ} - \tilde{D}_{ocean,poles})\sin(\tilde{\theta})^5& \text{for the ocean} \\
\tilde{D}_{NP} + (\tilde{D}_{equ} - \tilde{D}_{NP})\sin(\tilde{\theta})^5  &  \text{for the Northern Hemisphere} \\
\tilde{D}_{SP} + (\tilde{D}_{equ} - \tilde{D}_{SP})\sin(\tilde{\theta})^5  & \text{for the Southern Hemisphere} \\
\end{cases}
\end{align}

The corresponding thermal conductivity coefficient[^1] are as follows:

| **Types** | **Thermal conductivity in $W/K/m^2$** |
|-------------|-----------|
| Ocean at Poles | 0.4 |
| Ocean at Equator | 0.65 |
| Equator | 0.65 |
| North Pole | 0.28 |
| South Pole | 0.2 |

Note: For the surface type dependency of the heat conduction, you will once again need to include the function *read\_geography* from milestone 1.


3. Use the above formulas, the function *calc\_diffusion\_coefficients* and the functions for the other parameters from previous milestones to write a function *compute\_equilibrium\_EBM\_2D* to compute the EBM for two spatial dimensions using a forward Euler time integration scheme. Try to calculate the EBM using this function. What happens?

4. Write a function *calc\_jacobian\_EBM\_2D* that computes the Jacobian matrix of the FD operator numerically. You can calculate the $j$-th column of the Jacobian matrix with the equation:

\begin{equation*}
\underline{\mathbf{A}}_{\text{EBM}} \hat{\mathbf{K}}^j = \mathbf{R}(\hat{\mathbf{K}}^j)-\mathbf{R}(\mathbf{0})
\end{equation*}
where $\hat{\mathbf{K}}^j = (\hat{k}^j_1, \hat{k}^j_2, \ldots, \hat{k}^j_{NDOF})^T$ is a vector whose entries are defined as
\begin{equation}
\hat{k}^j_i =
\begin{cases}
1, & \mathrm{if} \  i = j, \\
0, & \mathrm{otherwise}.
\end{cases}
\end{equation}
 Note that the matrix $F$ does not depend on $T$ and cancels out.

5. Examine the sparsity pattern of the Jacobian matrix. What does the matrix look like?

6. Compute the eigenvalues of the Jacobian matrix. Determine the largest time step for which the scheme is stable using the forward Euler time integration scheme using the CFL condition.

[^1]: K. Zhuang, G.R. North, M.J. Stevens, _A NetCDF version of the two-dimensional energy balance model based on the full multigrid algorithm_, SoftwareX, Vol. 6, pp. 198-202, July 7, 2017.