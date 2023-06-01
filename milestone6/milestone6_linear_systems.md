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

**TODO: Add example**

@@colbox-blue
**Note:** We can convert a dense matrix to sparse format in Python using the SciPy 2-D sparse array package for numeric data:
```python
from scipy import sparse
sparse_jacoban = sparse.csc_matrix(dense_jacobian)
```
Similarly, we can convert a dense matrix to sparse format in Julia using the function `sparse` from the library `SparseArrays`:
```julia
using SparseArrays: sparse
sparse_jacoban = sparse(dense_jacobian)
```
@@

## Solving Sparse Linear Systems
There is a myriad of methods to solve linear systems efficiently
Solving with LU
