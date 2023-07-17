+++
title = "Milestone 1"
hascode = false
rss = "Description"
rss_title = "Milestone 1"
rss_pubdate = Date(2019, 5, 1)

tags = ["climatesystem"]
+++

# Milestone 1 - Geography and Visualization


## Read-in and Plot Geography

[(Download description in PDF format)](/assets/milestone1/description.pdf)\\
[(Download *The_World128x65.dat*)](/assets/milestone1/input/The_World128x65.dat)

In this milestone, the goal is to read-in the geography information of planet Earth and plot it as a two-dimensional world map.
We consider an equirectangular grid of the Earth, i.e., an equidistant rectangular grid in spherical coordinates, where the grid point $(i, j)$ has the spherical coordinates $(\theta_i, \varphi_j)$, where $\theta_i$ is the latitude between $-90\degree$ south and $90\degree$ north (including the poles) and $\varphi_j$ the longitude between $-180\degree$ west and $180\degree$ east.
The basis for this is the input file
[The_World128x65.dat](/assets/milestone1/input/The_World128x65.dat), which describes the distribution of the different Earth surface types.
This file contains a matrix $G \in \mathbb{N}^{65 \times 128}$ with entries $g_{ij} \in \{1,2,3,5\}$, where the entry $g_{ij}$ stores the Earth surface type at grid point $(i, j)$.
Here, $1$ represents the Earth surface type *land*, $2$ represents *sea ice*, $3$ represents *snow*, and $5$ represents *ocean*. The grid resolution in longitude and latitude direction is $2.8125^{\circ}$.
The basis for this distribution and grid is from Zhuang et al.[^1]
You can proceed as follows:

1. Write a function *read_geography*, which reads the file [The_World128x65.dat](/assets/milestone1/input/The_World128x65.dat) from the folder *input* and outputs a matrix $T \in \mathbb{N}^{65 \times 128}$  with the classification of the earth surface types.

2. Write a function *robinson_projection*, which maps an equirectangular grid in spherical coordinates to the plane.
For simplicity use the approximate formula by Beineke for the Robinson projection,
\begin{align*}
	&x\left(\theta, \varphi \right) = \frac{\varphi}{\pi}\left(0.0379\hspace{0.05cm}\theta^6 - 0.15\hspace{0.05cm}\theta^4 - 0.367\hspace{0.05cm}\theta^2 + 2.666 \right),\\
	&y\left( \theta, \varphi \right) = 0.96047\hspace{0.05cm}\theta - 0.00857\hspace{0.05cm}\text{sign}\left(\theta\right)\left|\theta\right|^{6.41},
\end{align*}
where $\varphi$ is the longitude and $\theta$ the latitude in radians.
This function should return two matrices $X = x_{ij}$ and $Y = y_{ij}$, where $x_{ij} = x(\theta_i, \varphi_j), y_{ij} = y(\theta_i, \varphi_j)$.

![](/assets/milestone1/Robinson_projection.jpg)
* The world on Robinson projection. Source: [Wikipedia](https://en.wikipedia.org/wiki/File:Robinson_projection_SW.jpg) ([Daniel R. Strebe](https://commons.wikimedia.org/wiki/User:Strebe)), 15th August 2011, [CC-BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/).

3. Write a function *plot_geo* that creates a plot of the Earth surface type $g_{ij}$ against the mapped coordinates $(x_{ij}, y_{ij})$.
4. Use these functions in a program and run it to check your results.

[^1]: K. Zhuang, G.R. North, M.J. Stevens, _A NetCDF version of the two-dimensional energy balance model based on the full multigrid algorithm_, SoftwareX, Vol. 6, pp. 198-202, July 7, 2017.




