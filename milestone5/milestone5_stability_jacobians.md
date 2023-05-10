+++
title = "Milestone 5"
hascode = true
rss = "Description"
rss_title = "Milestone 5"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 5 - Stability and Jacobian Matrices

\toc

In this chapter, we are interested in analyzing the stability properties of our [semi-discrete $2D$ EBM](/milestone5/milestone5_spatial_discretization/#semi-discrete_ebm).


## Mapping from a $2D$ Operator to a $1D$ Vector

As a first step, we will first rewrite our [system of ODEs in matrix form](/milestone5/milestone5_spatial_discretization/#eqsemidisc_ebm_matrix) as a system of ODEs in vector form, as that will facilitate the stability analysis of the system.

In particular, we want to concatenate the columns of our operators that are defined in matrix form (with $\nlat \times \nlong$), such that we create a one-dimensional vector with $\ndof = \nlat \times \nlong$ entries, where $\ndof$ stands for _number of degrees of freedom_ (see figure).

@@im-30
![](/assets/milestone5/mapping.png)
@@

Applying this mapping to our [system of ODEs in matrix form](/milestone5/milestone5_spatial_discretization/#eqsemidisc_ebm_matrix), we obtain
$$\label{eq:semidisc_ebm_vector}
\dot{\mathbf{T}} = 
\mathbf{R}(\mathbf{T})
+
\mathbf{F}(t),
$$
where $\mathbf{T}, \mathbf{R}(\mathbf{T}), \mathbf{F}(t) \in \R^{\ndof}$.

@@colbox-blue
**Remark:** To map a matrix ($2D$) variable to a $1D$ vector, you can use the command `vec` in Julia, and the command `.flatten()` in numpy (Python).
$$
\mat{R} \rightarrow \mathbf{R}
$$
@@

## Stability Theory for Linear Systems of Equations

As we saw in [Milestone 3](/milestone3/milestone3_odesolvers/#stability_of_the_explicit_and_implicit_euler_methods), the maximum allowable time-step size depends on the choice of the time discretization method for scalar ODEs.
In this section, we will apply the stability theory presented in [Milestone 3](/milestone3/milestone3_odesolvers/#stability_of_the_explicit_and_implicit_euler_methods) to systems of ODEs that result from spatial discretizations to study the stability properties of our $2D$ EBM model.

As a first step for our analysis, we will take advantage of the linearity of the term $\mathbf{R}(\mathbf{T})$ to rewrite our system of ODEs \eqref{eq:semidisc_ebm_vector} as
$$\label{eq:system_linear}
\dot{\mathbf{T}} = 
\underbrace{
    \mat{A} \mathbf{T}
}_{\mathbf{R}(\mathbf{T})}
+
\mathbf{F}(t),
$$
where $\mat{A}$ is the so-called Jacobian matrix of the system.

Note that \eqref{eq:system_linear} differs from the [linear scalar homogeneous ODE](/milestone3/milestone3_odesolvers/#eqsimpleode) that we analyzed in Milestone 3. As a matter of fact, \eqref{eq:system_linear} is nonhomegeneous and the different ODEs are coupled through matrix $\mat{A}$.

Assuming that our Jacobian matrix is invertible, we can rewrite it as 
$$
\mat{A} = \mat{V} \, \mat{\Lambda} \, \mat{V}^{-1},
$$
where the dense matrix $\mat{V} = [\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_{\ndof}]$ contains the eigenvectors of $\mat{A}$ in its columns, and $\mat{\Lambda} = \text{diag} (\lambda_1, \lambda_2, \ldots, \lambda_{\ndof})$ is a diagonal matrix that contains the eigenvalues of $\mat{A}$.
Hence, \eqref{eq:system_linear} is equivalent to the system of ODEs,
$$\label{eq:system_linear_decoupled}
\dot{\mathbf{Z}} = 
\mat{\Lambda} \mathbf{Z}
+
\pmb{\mathcal{F}}(t),
$$
where $\mathbf{Z} = \mat{V}^{-1}\mathbf{T}$ and $ \pmb{\mathcal{F}}(t) \coloneqq \mat{V}^{-1} \mathbf{F}(t)$.

Note that the individual variables in \eqref{eq:system_linear_decoupled} are decoupled due to the diagonal matrix $\mat{\Lambda}$. 
As a result, \eqref{eq:system_linear_decoupled} resembles the [linear scalar homogeneous ODE](/milestone3/milestone3_odesolvers/#eqsimpleode) and the solutions have the form
$$\label{eq:sol_diag}
Z_i(t) = a_i e^{\lambda_i t} + Z^p_{i}(t),
$$
where $Z^p_{i}$ is a so-called _particular solution_ (also _particular integral_) for the $i^{\text{th}}$ nonhomogeneous problem that depends on $\mathcal{F}_i(t)$ and $\lambda_i$, and $a_i$ is a constant that is adjusted to match the initial condition,
$$
a_i = Z_i(t=0) - Z^p_{i}.
$$

In the case that
$$
\lambda_i < 0, \qquad \text{for}~ i=1, 2, \ldots, \ndof,
$$
we expect that $\mathbf{Z}(t)$ converges to the _particular solution_ $\mathbf{Z}^p(t)$ as $t \rightarrow \infty$.

@@colbox-blue
**Remark:** The stability properties of our system of ODEs \eqref{eq:system_linear} are inherited from the linearized system \eqref{eq:system_linear_decoupled}.
As a matter of fact, we can compute the solutions of \eqref{eq:system_linear} from \eqref{eq:sol_diag} as
$$\label{eq:sol_linear}
\mathbf{T}(t) = \mat{V} \mathbf{Z}(t) = \sum_{i=1}^{\ndof} a_i e^{\lambda_i t} \mathbf{v}_i + \underbrace{\mat{V}\mathbf{Z}^p(t)}_{\mathbf{T}^{p}(t)},
$$
where $\mathbf{T}_{p}$ is now the _particular_ or _equilibrium_ solution of the original nonhomogeneous problem.
@@

We can now apply the stability theory presented in [Milestone 3](/milestone3/milestone3_odesolvers/#stability_of_the_explicit_and_implicit_euler_methods) to each scalar ODE of our diagonalized system \eqref{eq:system_linear_decoupled}.
For systems of ODEs, we require that **all** quantities $\lambda_i \Delta t$, for $i=1, 2, \ldots, \ndof$, lie in the stability region of the time-integration method used.

@@colbox-blue
**Remark:** The maximum allowable time-step size to solve a system of ODEs that results from a spatial discretization of a PDE depends on the spatial discretization method.
For instance, note that the eigenvalues $\lambda_i$ of matrix $\mat{A}$ depend on the spatial discretization method and grid size.
@@

## Computing Jacobian Matrices (Numerical Jacobian)

It is possible to expand a general nonlinear operator, $\mathcal{R}(\mathbf{T})$, with a multi-dimensional Taylor series to obtain,
$$ \label{eq:JacTaylorExpansion}
\mathcal{R}(\mathbf{T}) =\mathcal{R}(\mathbf{T}_0) + \mat{A(\mathbf{T}_0)} \Delta \mathbf{T}
+ \mathcal{O} ((\Delta \mathbf{T})^2),
$$
where the Jacobian matrix is defined as
$$
\mat{A} (\mathbf{T}_0) \coloneqq 
\left. \partialderiv{\mathcal{R} (\mathbf{T})}{\mathbf{T}} \right|_{\mathbf{T}_0}.
$$


This Taylor expansion can also be used to approximate the Jacobian matrix. For instance, one can assume that $\mathbf{T}$ differs from $\mathbf{T}_0$ only in one degree of freedom by a small quantity $\epsilon$, 
\begin{equation}
\mathbf{T} = \mathbf{T}_0 + \epsilon \hat{\mathbf{K}}^j,
\end{equation}
where $\hat{\mathbf{K}}^j = (\hat{k}^j_1, \hat{k}^j_2, \ldots, \hat{k}^j_{\ndof})^T$ is a vector whose entries are defined as
\begin{equation}
\hat{k}^j_i = 
\begin{cases} 
1 & \mathrm{if} \  i = j \\
0 & \mathrm{otherwise}.
\end{cases} 
\end{equation}

Equation \eqref{eq:JacTaylorExpansion} can then be rearranged as
\begin{equation} \label{eq:NumJacOrder1}
\mat{A} \hat{\mathbf{K}}^j = \frac{\mathcal{R}(\mathbf{T}_0+\epsilon \hat{\mathbf{K}}^j)-\mathcal{R}(\mathbf{T}_0)}{\epsilon} + \mathcal{O}(\epsilon ^ 2).
\end{equation}

For general nonlinear operators, the quantity $\mat{A} \hat{\mathbf{K}}^j$ is a good approximation for the $j^{\mathrm{th}}$ column of the Jacobian matrix if $\epsilon$ is small enough.

Since the right-hand-side term of our $2D$ EBM is a linear operator, \eqref{eq:NumJacOrder1} is exact for any $\epsilon$ and any linearization temperature $\mathbf{T}_0$, and the second and higher-order terms vanish. As a result, we can write a column of the Jacobian of our $2D$ EBM as
$$ \label{eq:NumJacEBM}
\mat{A}_{\text{EBM}} \hat{\mathbf{K}}^j = \mathbf{R}(\hat{\mathbf{K}}^j)-\mathbf{R}(\mathbf{0}).
$$

The whole Jacobian matrix can be recovered by computing \eqref{eq:NumJacEBM} $\ndof$ times (i.e., $\ndof + 1$ computations of the right-hand-side operator $\mathbf{R}(\mathbf{T})$):
$$
\mat{A}_{\text{EBM}} = \left[ \mat{A}_{\text{EBM}} \hat{\mathbf{K}}^{1}, \mat{A}_{\text{EBM}} \hat{\mathbf{K}}^{2}, \ldots, \mat{A}_{\text{EBM}} \hat{\mathbf{K}}^{\ndof}  \right].
$$


