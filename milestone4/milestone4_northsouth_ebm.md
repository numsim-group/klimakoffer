+++
title = "Milestone 4 - North/South EBM"
hascode = true
rss = "Description"
rss_title = "Milestone 4 - North/South EBM"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 4 - North/South EBM

\toc

As a next step, we want to investigate the difference of the southern and nothern hemisphere. If we visualize again the surface of Earth we
see a strong difference regarding the geography:

\fig{/assets/milestone4/map.png}
* Figure from Blue Marble Next Generation, NASA, 2004, Public Domain.

An obvious difference is that there is more ocean area in the south and more land area in the north. This is an important difference for our EBM
as the heat capacity coefficient $C(x)$ has strong variations: in particular, ocean grid cells (water) have a much larger value of $C$, which indicates that the thermodynamic time scales are much longer, which means that the ocean reacts much "slower" compared to the soil. 

It is even possible to estimate the thermal radiation relaxation timescales by inspecting again the temperature depending parts of our 0D-EBM,
$$
C\,T_t + A(CO_2) + B\,T = \text{solar forcing}, 
$$
where a simple physical units analysis shows that the ratio of $C$ and $B$ has the units of time, i.e., seconds $[s]$ (or years, depending on the temporal scaling of $C$)
$$
\tau := C / B. 
$$
The timescale $\tau$ is the relaxation time and is about $5.9$ years for the ocean and about $5.4$ days for the land. Hence, following these arguments, we do expect to see a difference in the annual temperature distribution between south and north. 

To investigate this behaviour we apply our time dependent EBM again, but not for the whole Earth. We are going to use one 0D-EBM for the southern hemisphere and one for the northern hemisphere. Hence, our goal is to get two OD-EBMs: for one we have to do area-averaging over the northern hemisphere
$$
\overline{C}_N\,\frac{\partial T_N}{\partial t}(t) + A(CO_2) + B\, T_N(t) = \overline{S}^{sol}_N(t),
$$
and for the other we have to do area-averaging over the southern hemisphere
$$
\overline{C}_S\,\frac{\partial T_S}{\partial t}(t) + A(CO_2) + B\, T_S(t) = \overline{S}^{sol}_S(t).
$$
To compute the special area averages, we remind ourselves about the computational grid and the definition of the area vector that we use to compute the area average:

\fig{/assets/milestone4/Grid.png}

We are now marking the equator with a blue dashed line.
We can see that we have to treat the grid points on the equator carefully, as their values contribute to the area averages of both the northern and southern hemispheres. 
If we use an odd number of points in the latitude direction, one half of the value of each equator grid cell contributes to the North and the other half to the South.

@@colbox-blue
**Remark:** We note that our grid always has an odd number of points in latitude direction. We can thus get the index of the equator grid cells by
\begin{align*}
j_{\text{equator}} &= (n_{latitude} - 1)/2 + 1 & \text{(1-based indexing)},
\\
\text{or} \qquad \qquad & &
\\
j_{\text{equator}} &= (n_{latitude} - 1)/2 & \text{(0-based indexing)}.
\end{align*}
@@

Let us assume that we have a 2D data field $F[j,i]$, with $j=1,...,\nlat$ and $i=1,...,\nlong$, for which we want to compute an area average over the nothern hemisphere. 
We assume that we have our area vector $area[j]$, with $j=1,...,\nlat$, from [milestone 3](/milestone3/milestone3_averages/). 
With the assumption that along the longitude direction $i$, all values at the north pole of $F[1,i]$ are the same, we can compute the area average as 
\begin{align}
\overline{F}_N = 2 \Bigg(& \texttt{area}[1]\,F[1,1]  & (\text{North~pole}) 
\\
+& \sum_{j=2}^{j_{\text{equator}}-1} \texttt{area[j]}\,\sum\limits_{i=1}^{\nlong} F[j,i]
& (\text{Intermediate~cells})
\\
&\left.+ \frac{1}{2}\,\texttt{area}[j_{\text{equator}}]\,\sum\limits_{i=1}^{\nlong} F[j_{\text{equator}},i] \right)
& (\text{Equator})
\end{align}
Note that we are multiplying by the factor $2$ as our $\texttt{area}[j]$ vector is normalized for the whole sphere surface, but the North is only half of it.
In an analogous way, we can easily compute the area average of the southern hemisphere.

With these new averages, and the tools we developed before, we are able to now generate a 0D-EBM for the northern and the southern hemispheres each. We solve the two 0D-EBM models individually to equilibrium as an approximation. The next figure shows the annual temperature distributsions: 


\fig{/assets/milestone4/temperatures_northsouth.png}

We can clearly see that the temperature distributions of north and south are phase shifted.  It seems, however, that the phase shift in this simple model is mainly due to the seasonal difference of the solar insolation. There is no clear effect due to the difference in heat capacity visible. The average heat-capacities are $\overline{C}_N = 5.2$ $[J/Ksm]$ and $\overline{C}_S=7$ $[J/Ksm]$. If we check the maximum temperature variation in the hemispheres, we get $\Delta T_{N,max} = 5$ and $\Delta T_{S,max} = 4.2$ $[K]$, which shows that a lower heat-capacity allows for a stronger variation in temperature. From our own experience about the annual seasons, it is clear that the seasonal temperature distribution from both, the North and the South, are still off. It seems that averaging our EBM first in the North and in the South, then compute the temperatures is too crude of an approximation to capture the seasonal temperature distribution accurately.
