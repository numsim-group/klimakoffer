+++
title = "Milestone 4"
hascode = true
rss = "Description"
rss_title = "Milestone 4"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 4

\toc

## Refining the 0D-EBM - Part I

As a next step, we want to investigate the difference of the southern and nothern hemisphere. If we visualize again the surface of Earth we
see a strong difference regarding the geography:

> INSERT MISSING STUFF HERE

An obvious difference is that there is more ocean area in the south and more land area in the north. This is an important difference for our EBM
as the heat capacity coefficient $C(x)$ has strong variations: in particular, ocean grid cells (water) have a much larger value of $C$, which indicates that the thermodynamic time scales are much longer, which means that the ocean reacts much "slower" compared to the soil. 

It is even possible to estimate the thermal radiation relaxation timescales by inspecting again the temperature depending parts of our 0D-EBM 
$$
C\,T_t + A(CO_2) + B\,T = \text{solar forcing}, 
$$
where a simple physical units analysis shows that the ratio of $C$ and $B$ has the units of time, i.e., seconds $[s]$
$$
\tau := C / B. 
$$
The timescale $\tau$ is the relaxation time and is about $4.3$ years for the ocean and about $31$ days for the land. Hence, following these arguments, we do expect to see a difference in the annual temperature distribution between south and north. 

To investigate this behaviour we apply our time dependent on EBM again but not for the whole Earth, but one 0D-EBM for each part, the southern hemisphere and one for the North. Hence, our goal is to get two OD-EBM: one area averaged over the northern hemisphere
$$
C_N\,\frac{\partial T_N}{\partial t}(t) + A(CO_2) + B\, T_N(t) = \overline{Q\alpha\,S}_N(t),
$$
and one area averaged over the southern hemisphere
$$
C_S\,\frac{\partial T_S}{\partial t}(t) + A(CO_2) + B\, T_S(t) = \overline{Q\alpha\,S}_S(t).
$$
To compute the special area averages, we remind ourselves about the computational grid and the definition of the area vector that we use to compute the area average:

> INSERT MISSING STUFF HERE

We can see that we have to be careful for the equator contribution where one half of the equator grid cells go to the North and the other half to the South.

@@colbox-blue
**Remark 1:** We note that our grid always has an odd number of points in latitude direction. We can thus get the index of the equator grid cells by
$$
j_{equator} = (n_{latitude} - 1)/2.
$$
@@

Let us assume that we have a 2D data field $F[j,i]$, with $j=0,...,n_{latitude}$ and $i=0,...,n_{longitude}$, for which we want to compute and area average, but only over the nothern hemisphere. We assume that we have our area vector $area[j]$, with $j=0,...,n_{latitude}$, from milestone 3. With the assumption that along the longitue direction $i$, all values at the north pole of $F[0,i]$ are the same, we can compute the first contribution of the north area average from the north pole as 
$$
\overline{F}_N = area[0]\,F[0,0]. 
$$
As discussed, the equator contribution only contributes with a half to the north average 
$$
\overline{F}_N =\overline{F}_N + \frac{1}{2}\,area[j_{equator}]\,\sum\limits_{i=0}^{n_{longitude}} F[j_{equator},i].
$$
The rest of the contributions are from inner grid nodes in the northern part:

> INSERT MISSING STUFF HERE

We need to be careful at the end, as our $area[j]$ vector is normalized for the whole sphere surface, but the North is only half of it 
$$
\overline{F}_N = \overline{F}_N / 2.
$$
In an analogous way, we can easily compute the area average of the southern hemisphere.

With these new averages, and the tools we developed before, we are able to now generate a 0D-EBM for the northern and the southern hemispheres each. We solve the two 0D-EBM models individually to equilibrium as an approximation. The next figure shows the annual temperature distributsions: 

> INSERT MISSING STUFF HERE
> **We should maybe include the temporal evolution of the source term S(t) in the above figure to investigate the phase shift regarding**
> I am not yet happy with the part that follows, to clearly explain, why the temperatures are not phase shifted...

We can clearly see that the temperature distributions of north and south are phase shifted.  It seems however, that the phase shift in this simple model is mainly due to the seasonal difference of the solar insoluation. There is no clear effect due to the difference in heat capacity visible. The average heat-capacities are $C_N = 5,2$ and $C_S=7$ $[J/Ksm]$. If we check the maximum temperature variation in the hemispheres, we get $\Delta T_{N,max} = 5$ and $\Delta T_{S,max} = 4.2$ $[K]$, which shows that a lower heat-capacity allows for a stronger variation in temperature. From our own experience about the annual seasons, it is clear that the seasonal temperature distribution from both, the North and the South, are still off. It seems that averaging our EBM first in the North and in the South, then compute the temperatures is to crude of an approximation to capture the seasonal temperature distribution accurately.

## Refining the 0D-EBM - Part II

As a next step, we want to actually use the local information $C(x)$, $\alpha(x)$, and $S(x,t)$ that we have available, see milestone 2. Hence, the idea for the next step improvement of our 0D-EBM is to set up local ODE problems separately for each grid node $j,i$. Then solve all the ODEs for each grid node to get a full 2D temperature field (3D field, if we count the time dependence as well). 

The local 0D-EBM for a local grid node $i...longitude$ and $j...latitude$ is given by 
$$
C_{ji}\,\frac{\partial T_{ji}}{\partial t}(t) + A(CO_2) + B\, T_{ji}(t) = Q\,(1-\alpha_{ji})\,S_{ji}(t),
$$
and is solved separately for each node. After solving the ODEs, we can then average the resulting time dependend temperature field $T[j,i,k]$, over the whole globe, and/or separately for the norther and the southern hemispheres. The results are plotted in the next figure

> INSERT MISSING STUFF HERE

We can observe some significant changes in the quality of our approximation: 

(i) The difference in the temperature variations in the North and South is much bigger now with $\Delta T_{N,max} \approx 33$ and $\Delta T_{S,max} \approx 17$ $[K]$.

(ii) Both, the North and South temperature distributions are now phase shifted compared to our 0D-EBM North-South from before. Compared to the real observed temperature distribution (and our own experience), these are much more realistic. We now can even see that the South is a bit "later" than the North, which fits our expectation because of the difference of ocean area in the South and North. 

We encourage the interested reader to try out different values of $CO_2$ and to compare the data to real measured temperatures. In particular, the sensitivity of changing the $CO_2$ value and its effect on the global average temperature change $\Delta t$.

It is clear that the local 0D-EBM for each grid node has in general improved the solution quality and its physical realism. However, it is still far off from observations. In particular, the minimum and maximum temperatures seem a bit too extreme. Furthermore, the solution field appears discontinuous across strong geography changes, which is also not realistic. 

As an intermediate conclusion, it seems that a major process is missing: a process that should smear out the extremes a bit (reducing minimum and maximum values) and that should also connect the different grid node solutions with each other, such that neighbour node temperatures are somewhat continuous and to get rid of the very steep gradients across geography changes. 

Such a process, such a behaviour is precisely the effect/role of diffusion and dissipation. Hence, the next step to improve our 0D-EBM approach is to move away and beyond the ODE modeling and consider a PDE by incorporating heat transfer and the heat flux, i.e., the diffusion operator $\vec{\nabla}\cdot(D(x)\,\vec{\nabla} T)$.
