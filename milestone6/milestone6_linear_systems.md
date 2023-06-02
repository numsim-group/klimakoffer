+++
title = "Milestone 6"
hascode = true
rss = "Description"
rss_title = "Milestone 6"
rss_pubdate = Date(2022, 5, 1)

tags = ["linear systems", "linear solvers"]
+++

# Milestone 6 - Solving Systems of Linear Equations

\toc

## Sparse Matrix Format

The spatial discretization that we derived in [milestone 5](/milestone5/milestone5_spatial_discretization) is compact, which means that the operators only couple the solution at each degree of freedom $j,i$ with the solutions at neighboring degrees of freedom: $L_{j,i} = L_{j,i} (T_{j,i}, T_{j,i+1}, T_{j,i-1}, T_{j+1,i}, T_{j-1,i})$.
As a consequence, we obtain a very sparse system matrix $\mat{M}$ (many of its entries are zeros).
The sparsity of the system matrix is exactly the same as the sparsity of $\mat{A}$ that we computed in [the assignment of milestone 5](/milestone5/milestone5_results/#sparsity_pattern_of_the_jacobian) because we are only adding values to the diagonal entries of $\mat{A}$:
$$
\mat{M} = \mat{I} - \Delta t \, \mat{A}.
$$

Utilizing a sparse matrix format can significantly enhance efficiency due to the abundance of zero entries in $\mat{M}$. This format optimizes storage space allocation by only storing the non-zero elements.

Sparse matrix formats offer additional advantages beyond storage savings. They can accelerate various matrix operations. To illustrate this, consider a matrix-vector multiplication involving a dense matrix and a dense vector. In this case, $\ndof^2$ multiplications and $\ndof (\ndof - 1)$ additions are required, where $\ndof \times \ndof$ represents the matrix size.
In contrast, employing a sparse matrix and a dense vector for the matrix-vector multiplication requires only $nnz$ multiplications and $\mathcal{O}(nnz)$ additions, where $nnz$ denotes the number of non-zero elements.

@@colbox-blue
**Remark:** Our spatial discretization operator couples each interior node with itself and four other nodes, each node at the north pole with all the nodes at the north pole and all the nodes at a colatitude $\colat = \Delta \colat$, and each node at the south pole with all the nodes at the south pole and all the nodes at a colatitude $\colat = \pi - \Delta \colat$.

Therefore, the number of non-zeros of the matrix is 
\begin{align}
nnz &= 5 (\nlat - 2) \, \nlong + 4 \, \nlong^2
\\
&=
6.5 \, \nlong^2 - \nlong,
\end{align}
which is much smaller than
\begin{align}
\ndof^2 &= \nlong^2 \, \nlat^2
\\
&= \frac{1}{4} \nlong^4 + \nlong^3 + \nlong^2.
\end{align}
@@

There are many different sparse matrix formats available in the literature. 
Depending on the application, some formats can offer more advantages than others. 
In this course, we will use the _Compressed Sparse Row_ (CSR) format, which is also known as the _Compressed Row Storage_ (CRS) format.
In the CSR format we store a sparse matrix using only three one-dimensional vectors:
* $\texttt{V}$: A vector of real numbers of size $nnz$ that contains the values of the entries of $\mat{M}$ in ascending row and ascending column order.
* $\texttt{COL\_INDEX}$: A vector of integers of size $nnz$ that contains the column index for each entry in $\texttt{V}$.
* $\texttt{ROW\_INDEX}$: A vector of integers of size $\ndof + 1$ that contains the index in $\texttt{V}$ and $\texttt{COL\_INDEX}$ where the given row starts. The last entry of $\texttt{ROW\_INDEX}$ contains the fictitious index in $\texttt{V}$ after the last valid index. In other words, it contains $nnz$ in zero-based programming languages (e.g., Python) and $nnz+1$ in one-based programming languages (e.g., Julia).

@@colbox-blue
**Example:** We want to store the matrix $\mat{M}$ in sparse format:
$$
\mat{M} =
\begin{pmatrix}
1 & 0 & 0 & 2 & 0 \\
3 & 4 & 0 & 5 & 0 \\
6 & 0 & 7 & 8 & 9 \\
0 & 0 & 10 & 11 & 0 \\
0 & 0 & 0 & 0 & 12 
\end{pmatrix}
\in \R^{5 \times 5},
$$
where $\ndof = 5$ and $nnz = 12$. We have:
* In a zero-based language (e.g., Python):
\begin{align}
\texttt{V} &= \left( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 \right)
\\
\texttt{COL\_INDEX} &= \left( 0, 3, 0, 1, 3, 0, 2, 3, 4, 2, 3, 4 \right)
\\
\texttt{ROW\_INDEX} &= \left( 0, 2, 5, 9, 11, 12 \right).
\end{align}
* In a one-based language (e.g., Julia):
\begin{align}
\texttt{V} &= \left( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 \right)
\\
\texttt{COL\_INDEX} &= \left( 1, 4, 1, 2, 4, 1, 3, 4, 5, 3, 4, 5 \right)
\\
\texttt{ROW\_INDEX} &= \left( 1, 3, 6, 10, 12, 13 \right).
\end{align}
@@

@@colbox-blue
**Note:** We can convert a dense matrix to sparse format in Python using the SciPy 2-D sparse array package for numeric data:
```python
from scipy import sparse
sparse_jacoban = sparse.csr_matrix(dense_jacobian)
```
Similarly, we can convert a dense matrix to sparse format in Julia using the function `sparse` from the library `SparseArrays`:
```julia
using SparseArrays: sparse
sparse_jacoban = sparse(dense_jacobian)
```
@@

## Solving Sparse Linear Systems

There is a myriad of methods to solve linear systems efficiently. 

Direct methods are the most reliable methods for solving _small_ linear systems since they obtain an exact solution to the problem to within rounding/conditioning errors. 
The most well known direct method is Gaussian elimination, which reduces the system to a triangular one using row operations. 
Once the system is triangular, the solution can be found via backward substitution. 
Gaussian elimination is considered to be a very expensive method since it requires $\mathcal{O}(2 \ndof^3/3)$ floating point operations. 
Many other alternatives are available in the literature, such as the Cholesky decomposition (for Hermitian, positive-definite matrices), the QR decomposition and the LU decomposition (both for general nonsymmetric systems). In next section, LU decomposition is briefly described, as it is the most used and efficient method for real nonsymmetric systems.

### LU Decomposition

The idea of LU decomposition is to factorize the system matrix as the product of non-singular lower- and upper-diagonal matrices $L$ and $U$,
\begin{equation}
\mat{M} = \mat{L} \, \mat{U}.
\end{equation}

Depending on the structure of matriz $\mat{M}$, this is not always possible. However, it a proper reordering of the matrix by rows or columns, also called permutation, is sufficient to enable the LU factorization of general nonsingular matrices,
\begin{equation}
\mat{P} \mat{M} = \mat{L} \, \mat{U},
\end{equation}
where $\mat{P}$ is known as the permutation matrix.

Once the matrix is factorized, the linear system $\mat{M} \mathbf{T}^{n+1} = \mathbf{Y}^{n+1}$ can be solved in two steps:
1. Solve the linear system $\mat{L} \mathbf{Z} = \mat{P} \mathbf{Y}^{n+1}$ for $\mathbf{Z}$.
2. Solve the linear system $\mat{U} \mathbf{T}^{n+1} = \mathbf{Z}$ for $\mathbf{T}^{n+1}$.

Since the matrices $\mat{L}$ and $\mat{U}$ are triangular, both linear solves can be done by backward and forward substitution. Therefore, the solving process is very cheap, requiring approximately $\mathcal{O}(\ndof^2/2)$ operations. The expensive part is the factorization of the matrix, which requires around $\mathcal{O} (\ndof^3/3)$ operations. 
As a consequence, the LU factorization is a feasible possibility when the factorized matrix can be used for several linear solves with different right-hand-sides, and there is enough storage for $\mat{L}$, $\mat{U}$, $\mat{M}$ and the variables needed to factorize the original system.

@@colbox-blue
**Note:** The LU decomposition of a sparse matrix can be done in Julia using the `LinearAlgebra` library. To perform the factorization, use the function `lu`:
```julia
using LinearAlgebra: lu
lu_decomposition = lu(sparse_jacobian)
```
After this command, the variable `lu_decomposition` stores the factorized matrices $L$ and $U$.
To solve the system for a particular right-hand side `rhs` and store the solution in the array `sol`, use the function `ldiv!`:
```julia
using LinearAlgebra: ldiv!
ldiv!(sol, lu_decomposition, rhs)
```

Similarly, LU decomposition of a sparse matrix can be done in Python using the subpackage `linalg` of `scipy.sparse`. 
To perform the factorization, use the function `factorized`:
```python
solve = sparse.linalg.factorized(sparse_jacobian)
```
To solve the system for a particular right-hand side `rhs` and store the solution in the array `sol`, we can use the newly defined function `solve`:
```python
sol = solve(rhs)
```
@@