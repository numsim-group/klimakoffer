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

As a first step, will first rewrite our [system of ODEs in matrix form](/milestone5/milestone5_spatial_discretization/#eqsemidisc_ebm_matrix) as a system of ODEs in vector form, as that will facilitate the stability analysis for our system.

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

Extension of [stability theory presented in Milestone 3](/milestone3/milestone3_odesolvers/#stability_of_the_explicit_and_implicit_euler_methods)

Since our two-dimensional EBM model is linear, we can write our semi-discretization as
$$\label{eq:system_linear}
\dot{\mathbf{T}} = 
\underbrace{
    \mat{A} \mathbf{T}
}_{\mathbf{R}(\mathbf{T})}
+
\mathbf{F}(t),
$$
where $\mat{A}$ is the so-called Jacobian matrix of the system.

To analyze the stability properties of our model, it will be useful to rewrite \eqref{eq:system_linear} as a collection of [simple ODE](/milestone3/milestone3_odesolvers/#eqsimpleode). 

Assuming that our Jacobian matrix is invertible, we can rewrite it as 
$$
\mat{A} = \mat{Q} \mat{\Lambda} \mat{Q}^{-1}
$$

$$\label{eq:system_linear}
\mat{Q}^{-1} \dot{\mathbf{T}} = 
\mat{\Lambda} \mat{Q}^{-1} \mathbf{T}
+
\mat{Q}^{-1} \mathbf{F}(t).
$$

## Computing Jacobian Matrices

### Numerical Jacobian

It is possible to expand a general nonlinear operator, $\mathcal{R}(\mathbf{T})$, with a Taylor series to obtain,
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
\mat{A} \hat{\mathbf{K}}^j = \frac{\mathcal{R}(\mathbf{T}_0+\epsilon \hat{\mathbf{K}}^j)-\mathcal{R}(\mathbf{T}_0)}{\epsilon} + \mathcal{O}(\epsilon ^ 2).
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

### Bonus: Analytical Jacobian

