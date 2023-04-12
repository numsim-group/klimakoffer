+++
title = "Milestone 2"
hascode = true
rss = "Description"
rss_title = "Milestone 2"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 2 - Albedo, Solar Forcing and Visualization in Time

## Compute Albedo and Solar Forcing Terms and Visualize in Time

[(Download description in PDF format)](/assets/milestone2/description.pdf)

In milestone 2 you will program the calculation of all necessary parameters for the calculation of the time and space dependent solar forcing and then visualize it as an animation over a period of one year in 48 time steps. For this you will need the functions *read_geography* and *robinson_projection* from milestone 1.

1. Write a function *calc_albedo* which outputs the albedo $\alpha_{ij}$ for the points $\left(\lat_i,\varphi_j\right)$ as a matrix $A \in \mathbb{R}^{65\times 128}$ depending on the earth surface type $T\in \mathbb{N}^{65\times 128}$ (see milestone 1). For the calculation use the albedo scheme below as proposed by Zhuang et al. [^1].

| **Land mask** | **Albedo** |                                                                                                 
|---------------|--------------|
| Land | $0.30+0.12P_2(x)$ |
| Ocean | $0.29+0.12P_2(x)$ |                                                                                               
| Sea ice | $0.60$ |                                                                                               
| Snow | $0.75$ |                                                                                                     

For $P_2(x)=\frac{(3x^2-1)}{2}$ and $x=\sin(\lat)$, where $\lat$ is the latitude in radians.

2. In the table below, besides atmospheric parameters, constants are given for each type of the earth's surface, which are used for the calculation of the effective heat capacity. The effective heat capacity of a grid cell results from the sum of the heat capacity of the respective surface type and the heat capacity of the overlying atmospheric column. Since the climate model measures time in years, this heat capactiy has to be divided by the number of seconds in a year $3.15576\times10^7 s$. Implement a function *calc_heatcapacity* analogous to *calc_albedo* to calculate the geography-dependent heat capacity per surface area using the given constants from the table below.

| **Types** | **Density $\rho$ in $\mathrm{kg\,m}^{-3}$}** | **Specific heat capacity in $\mathrm{J\,kg^{-1}K^{-1}}$** | **Depth in $\mathrm{m}$** |
|-------------|-----------|---------------|---------------|
| Atmosphere | 1.225 | 1000 | 3850 |
| Ocean mixed layer | 1030 | 4000 | 70 |
| Sea ice | 917 | 2000 | 1.5 |
| Land | 1350 | 750 | 1.0 |
| Snow | 400 | 880 | 0.5  |

3. Implement a function *read\_true\_longitude* to read in the true longitude $\lambda$ for each time step from the input file [True_Longitude.dat](/assets/milestone2/input/True_Longitude.dat) and output it as a vector $\lambda \in \mathbb{R}^{48}$.

4. Write a function *insolation* that calculates the insolation $S\left(\lat,t\right)$ as in Berger[^2]. Input arguments of the function should be:

* $\theta$: latitude,
* $S_0$: solar constant,
* $e$: eccentricity,
* $\epsilon$: obliquity,
* $\tilde \omega$: precession distance,
* $\lambda$: true longitude.

The insolation at latitude $\lat$ and time $t$ is given by:
\begin{align*}
S(\lat,t) =
\begin{cases}
0 & \text{for $(\lat,t)$ without sunrise}, \\
S_0\rho(t) \sin (\lat) \sin \big(\delta(t)\big) & \text{for $(\lat,t)$ without sunset}, \\
\frac{S_0 \rho(t)}{\pi}\Big(H_0(t) \sin (\lat) \sin \big(\delta(t)\big) + \cos (\lat) \cos \big(\delta(t)\big) \sin \big(H_0(t)\big)\Big) & \text{else}, \\
\end{cases}
\end{align*}
where $\delta(t)$ is the declination of the sun with
\begin{align*}
\sin(\delta(t)) = \sin(\epsilon)\sin (\lambda(t)).
\end{align*}
Here, $\rho$ is the normalized distance between earth and sun squared,
 \begin{align*}
 \rho(t) =  \left(\frac{1 - e \cos (\nu(t))}{1-e^2}\right)^2,
 \end{align*}
with eccentricity $e$ and the positional angle of the Earth on its orbit around the Sun $\nu(t)=\lambda(t)-\tilde \omega$ as given counterclockwise from aphelion.
\begin{align*}
H_0(t)= \arccos\big(-\tan(\lat) \tan(\delta(t))\big)
\end{align*}
describes the absolute angle of the sun at sunrise and sunset.
Make sure that your implementation handles values $\pm \infty$ of the tangent at the poles correctly.

5. Write a function *calc\_solar\_forcing* with input arguments

* $\alpha$: albedo values at the grid points,
* $S_0$: solar constant with default value $S_0 = 1371.685$,
* $e$: eccentricity with default value $e = 0.01674$,
* $\epsilon$: obliquity with default value $\epsilon  = 0.409253$,
* $\tilde \omega$: precession distance with default value $\tilde \omega = 1.783037$,
* $\left(\lambda_j\right)_{j = 1}^{n_t} \in \mathbb{R}^{n_t}$: true longitude for each time step.


Note that these default orbital parameters correspond to the year 1950.
This function should calculate and return the solar forcing $S_{sol}(\lat_i,\varphi_j,t_k)$ at each grid point $(\lat_i, \varphi_j)$ and for each time step $t_k$.
The return value should be a three-dimensional array, where the first dimension corresponds to the latitude, the second to the longitude, and the third to the time step.
The solar forcing is given by
\begin{equation*}
	S_{sol}\left(\lat,\varphi,t\right) = S\left(\lat,t\right)\alpha_c\left(\lat,\varphi\right).
\end{equation*}
Here, $S$ is the insolation defined above and $\alpha_c$ is the coalbedo given by $\alpha_c = 1 - \alpha$, where $\alpha$ is the albedo.


6. To visualize the solar forcing as an animation, calculate the albedo and the solar forcing for a grid size of $65 \times 128$ (with the [data file](/assets/milestone2/input/The_World128x65.dat) from milestone 1) and $48$ time steps within a year. Then plot your result using the function
*robinson_projection* from milestone 1.

[^1]: K. Zhuang, G.R. North, M.J. Stevens, {\em A NetCDF version of the two-dimensional energy balance model based on the full multigrid algorithm}, SoftwareX, Vol. 6, pp. 198-202, July 7, 2017.
[^2]: A. L. Berger, {\em Long-Term Variations of Daily Insolation and Quaternary Climatic Changes}, Journal of the Atmospheric Sciences, Vol. 35, pp. 2362-2367, 1987.}




