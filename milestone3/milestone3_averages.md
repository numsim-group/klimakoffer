+++
title = "Milestone 3"
hascode = true
rss = "Description"
rss_title = "Milestone 3"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 3 - Averages

\toc

---

## Area average

We have to be careful about computing the
area average: while we use a Cartesian computational grid in longitude and latitude, we need the average over the
sphere surface.

We consider the latitude/longitude grid:

\fig{/assets/milestone3/Grid.png}

In the figure, the discretization nodes are marked in purple, the grid _lines_ are marked in gray, and we have drawn additional auxiliary grid lines (dotted lines in red) between our computational nodes.
We are going to use these auxiliary grid lines to define the surface area that corresponds to each computational node.

Observe that, even though we use a regular grid, the area that corresponds to each computational node of our grid depends only on the latitude, and not on the longitude.
Therefore, our goal is to compute the corresponding
approximation of the area for each latitude
grid coordinate ($\lat_j$).

In addition, observe that each node at the pole maps to a triangle, and not to a quadrilateral.
In preparation for the discretization of the diffusion operator in sperical
coordinates ([milestone 5](/milestone5/milestone5_menu/)), we already start to
treat the pole regions ($j=1$ and $j=\nlat$)
in a special way.

### Poles (North and South)

We need to compute the area of a spherical cap at the angle $\lat = \frac{\Delta \lat}{2}$.

\fig{/assets/milestone3/SphericalCap.png}
* Source: Wikipedia

Using mathematical analysis, we can derive the area of a spherical cap with radius $R_E$ as
$$
\tilde A_{\text{pole}} = 2 \pi R_E^2 (1-\cos \lat) = 2 \pi R_E^2 \left(1-\cos \left(\frac{\Delta \lat}{2}\right)\right).
$$

As we will use the area calculation only in to compute area averages of quantities on the surface of the spghere, we will directly normalize this value with the surface area of the Earth (spehere), i.e.
$$\label{eq:cap}
A_{\text{pole}} = \frac{1}{4\pi R_E^2} \tilde A_{\text{pole}} = \frac{1}{2} \left(1-\cos \left(\frac{\Delta \lat}{2}\right)\right).
$$

### Interior nodes

For the interior nodes, we consider a spherical segment between the two angles $\lat_j - \Delta \lat/2$ and $\lat_j + \Delta \lat/2$.
We can use \eqref{eq:cap} to compute the difference between two spherical caps and obtain
\begin{align}
A_{\text{int}} &= \frac{1}{2} \left(\cos \left(\lat_j - \frac{\Delta \lat}{2}\right) - \cos \left(\lat_j + \frac{\Delta \lat}{2}\right)\right)
\\
&= \frac{1}{2}  \sin \left(\frac{\Delta \lat}{2}\right) \sin \left(\lat_j\right).
\end{align}

### Implementation and data structures

To compute the area average of a vector containing values at the grid points, we will directly store the normalized area entroes in a vector of floats with length $\nlat$ that we will call `area`.

For instance, assuming $\nlat = 65$, we can declare the vector in julia as
```julia:./define_area.jl
n_latitude = 65
n_longitude = 2 * (n_latitude - 1)
area = zeros(Float64,n_latitude)
```

The grid size is 
$$
\Delta \lat = \frac{\pi}{\nlat - 1}.
$$
Therefore, taking into account that julia uses one-based indexing we have for the poles
```julia:./poles.jl
delta_theta = pi / (n_latitude - 1)
area[1] = 0.5 * (1 - cos(0.5 * delta_theta))
area[n_latitude] = area[1]
```

Next, we store the area of the interior grid points, however, _with a twist_:
```julia:./fillarea.jl
for j=2:n_latitude-1
    theta_j = delta_theta * (j - 1)
    area[j] = sin(0.5 * delta_theta) * sin(theta_j) / n_longitude
end
# We print the area array to check that everything is right
println(area)
```
To obtain:
\output{./fillarea.jl}

@@colbox-blue
**Remark:** In the last step, we further divided the area of the inner sphere segments by the number of (regular) longitude grid cells. This means that, while the first and last entries of the array `area` contain the area of the entire spherical cap, the other entries of `area` only contain the cell-wise areas.
The reason is that all the grid points at the pole are co-located. As a result, all entries of the first and last row must be the same for any matrix that stores grid data in our model:
$$
T(j=1,i=1) = T(j=1,i=2) = \cdots = T(j=1,i=\nlong)
$$
\begin{align}
T(j=\nlat,i=1) &= T(j=\nlat,i=2) 
\\
&= \cdots = T(j=\nlat,i=\nlong)
\end{align}
Hence, it makes sense to not differentiate between those values, but just use the first one.

On the other hand, the solution values in the interior of the grid are in general different for each location `i`, e.g.:
$$
T(j,i=1) \ne T(j,i=2) = \cdots = T(j,i=\nlong),
$$
for $j\notin\{1,\nlat\}$. 
That is why we store the areas of the individual cells (normalization with the factor $1/\nlong$) for the interior cells.
@@

@@colbox-blue
**Remark:** We have to adjust the expressions if we use a zero-based programming language (like python).
@@

With this, we get an easy formula for the calculation of the area average of a given field `F(j,i)` wirh $j=1, \ldots, \nlat$ and $i=1, \ldots, \nlong$.
```julia:./computearea.jl
function calc_mean(field, area, n_latitude, n_longitude)
  if size(field)[1] !== n_latitude || size(field)[2] !== n_longitude
    error("field and area sizes do not match")
  end

  # Initialize mean with the values at the poles
  mean = area[1]*field[1,1] + area[end]*field[end,end]

  for j in 2:n_latitude-1
      for i in 1:n_longitude
          mean += area[j] * field[j,i]
      end 
  end

  return mean
end 
```

We test our routine by computing the area average of the dara $F(j,i)=1 \forall (i,j)$ to get the value of $1$.
```julia:./testarea.jl
field = ones(Float64,n_latitude,n_longitude)
println("Area : ", calc_mean(field, area, n_latitude, n_longitude))
println("Error: ", calc_mean(field, area, n_latitude, n_longitude)-1)
```
When we execute the code, we obtain:
\output{./testarea.jl}

## Temporal average

The temporal average is less involved and relatively straightforward. We make the assumption that tht temporal behavior is periodic (Earth moving around the sun in one year).

We normalize the time, such that $t=1$ corresponds to one year time.

In our temporal discretization, we will choose a fixed number of time steps per year $\ntime$ and compute the fixed time-step size as $\Delta t = 1 / \ntime$.

Given a vector that contains solution values in time $F(k)$ with $k=1, \ldots, \ntime$, we compute the temporal average as
\begin{align}
\widehat{F} &= \frac{\left( \sum_{k=1}^{\ntime} F(k) \Delta t_k \right)}{\left( \sum_{k=1}^{\ntime} \Delta t_k \right)}
\\
&= \Delta t \sum_{k=1}^{\ntime} F(k).
\end{align}