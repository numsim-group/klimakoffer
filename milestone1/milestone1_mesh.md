+++
title = "Milestone 1"
hascode = false
rss = "Description"
rss_title = "Milestone 1"
rss_pubdate = Date(2019, 5, 1)

tags = ["climatesystem"]
+++

# Milestone 1 - Meshing the sphere
\toc

In numerical methods, we often partition the domain in smaller subdomains using what is called a mesh. The subdomains are often called grid cells, or mesh cells, or grid boxes. We will then approximate the solution to the partial differential equation that describes our model within each of those little grid cells.

The easiest way to partition a sphere is to use grid cells with uniform radius, latitude/colatitude and longitude spacings:

![](/assets/milestone1/PolesProblem.png)
* Spherical coordinate system. Source: [Wikipedia](https://commons.wikimedia.org/wiki/File:Spherical_coordinate_system.svg) ([User: SharkD](https://commons.wikimedia.org/wiki/User:SharkD)), [CC-BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

We can observe from the figure that the
grid cells get smaller the closer they are to the
poles. So a regular grid in (co-)latitude and
longitude space gives an irregular grid on the
sphere surface.
In fact, we can note that [the determinant](/milestone1/milestone1_sphere/#eqdet) can
get equal to zero for $\colat \in \{ 0, \pi\}$. From linear
algebra, we know that matrices with determinant
equal to zero are not regular, i.e., they cannot
be inverted. Transformations where the Jacobian
matrix gets irregular are singular at this specific
locations. The locations $\colat=0$ (North) and $\colat=\pi$ (South)
correspond to the poles. Therefore, the transformation is
singular at the poles, which shows that, regarding
mappings, the poles are special locations and
are somewhat problematic when trying to
mesh the surface.


Because of the issues discussed above,
there are many alternative ways of constructing
meshes for sphere:

![](/assets/milestone1/GridsSphere.png)
* **Source**:  Example of horizontal grid types used by atmospheric models. [Weather forecasting models](https://www.encyclopedie-environnement.org/en/air-en/weather-forecasting-models/), Sylvie Malardel, Encyclopedia of the Environment, [CC-BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).



For instance, the famous ICON (Icosahedral Nonhydrostatic) model from the German weather service (DWD: Deutscher Wetterdienst) uses triangle
surface grids:


![](/assets/milestone1/ICONgrid.png)
* **Source**: [www.dwd.de](http://www.dwd.de)

Although we have identified some issues with grids that are regular in latitude/colatitude and longitude, we have decided to use this type of grid in our model. Despite its limitations, it is the simplest grid to construct and work with. To address the singularities at the poles, we will develop and implement a special fix.

The grid that we use in this course is illustrated below:

![](/assets/milestone1/OurGrid.png)

In the illustration, the boundaries of our domain are marked in blue, the grid lines of the mesh are marked in gray and the grid points of the mesh (the positions where we will store our numerical solution) are marked as purple circles.

Since our domain is periodic in the longitude direction, we do not need to store the last column of grid points ($\varphi = \pi$), as their position on the surface of the sphere is the same as for the first column ($\varphi = -\pi$).

Note that all the points in the first row ($\lat = \pi/2$) correspond to the same position (north pole), and all the points in the last row ($\lat = -\pi/2$) too (south pole).
Nevertheless, we will keep these duplicated grid points in our model to simplify the storage (we can use a matrix format to store the values).

We define the number of grid points as
\begin{align}
\nlong \in \mathbb{N} \\
\nlat \in \mathbb{N}
\end{align}
to get the size of the grid cells as
\begin{align}
\Delta \varphi &= \frac{2\pi}{\nlong} \\
\Delta \lat &= \Delta \colat = \frac{\pi}{\nlat -1}.
\end{align}

To simplify the derivation of our numerical climate model, we will use uniform grids, i.e., $\Delta \varphi = \Delta \lat$, so we require:
$$
\nlong = 2(\nlat - 1).
$$

The longitude grid node locations are defined from west to east
\begin{align}
\varphi_i &= -\pi + (i-1) \Delta \varphi, & i&=1, 2, \ldots, \nlong,
\end{align}
the latitude grid node locations are defined from north to south
\begin{align}
\lat_j &= \frac{\pi}{2} - (j-1) \Delta \lat, & j&= 1, 2, \ldots, \nlat,
\end{align}
and equivalently the colatitude grid node locations are
\begin{align}
\colat_j &= (j-1) \Delta \colat, & j&= 1, 2, \ldots, \nlat.
\end{align}

To illustrate the results of the model,
we could directly plot in computational space (lat/long grid).
However, it is more pleasing to the eyes and
more common to plot the results in physical
space (or an approximation to that).
There are many options available in the literature
to get from lat/long to other coordinates.
In this course we use the so-called Robinson
projection, which is an interesting one as this
transform has no mathematical properties (it does not keep distances, angles, or areas) but
was designed by Robinson by hand to look
pleasing to his eyes(!).

While Robinson provided a translation
table for values (lat/long), there are
closed form approximations available.
The form we consider was presented by
> Beineke, D. (1991). Untersuchung zur Robinson-Abbildung und Vorschlag einer analytischen Abbildungsvorschrift. Kartographische Nachrichten, 41(3), 85-94.

The approximation that Beineke proposes
is a simple-to-evaluate algebraic formula.
Given the computational coordinates, he defines the new coordinates of the
curvilinear grid as
\begin{align}
    \hat{x}_{ij} &= (d + e \lat_j^2 + f \lat_j^4 + g \lat_j^6) \varphi_i \\
    x &= \frac{180}{\pi} \frac{\hat{x}_{ij}}{\max_{i,j} (\hat{x}_{ij})} \\
    \hat{y}_{ij} &= a \lat_j + b \texttt{sign} (\lat_j) |\lat_j|^c \\
    y &= 90 \frac{\hat{y}}{\max_{i,j} (\hat{x}_{ij})}
\end{align}
with the parameters
\begin{align}
a =&  0.96047, &
    b =& -0.00857, &
    c =&  6.41, \\
    d =&  2.6666, &
    e =& -0.367, &
    f =& -0.150 \\
    g =&  0.0379.
\end{align}

Applying the simplified Robinson
projection to the computational grid
gives a curvilinear mesh a shown in the
figure:

![](/assets/milestone1/Robinson.png)
* Map of the world in a Robinson projection with Tissot's Indicatrix of deformation. Source: [Wikipedia](https://commons.wikimedia.org/wiki/File:Tissot_indicatrix_world_map_Robinson_proj.svg) (Eric Gaba - [User:Sting](https://commons.wikimedia.org/wiki/User:Sting)), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

As can be seen in the figure, the
geographic map can be used to
display the land-sea-ice-show
mask of Earth, which is the
topic of the first assignment, milestone 1.

