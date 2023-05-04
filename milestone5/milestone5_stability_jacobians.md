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

## Stability Theory for Linear Systems of Equations

* TODO: Explain mapping from $\nlat \times \nlong$ operator to a $1D$ vector. Julia command `vec` and numpy (Python) commands `.flatten()` and `reshape`.
$$
\mat{R} \rightarrow \mathbf{R}
$$

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

