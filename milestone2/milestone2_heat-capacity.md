+++
title = "Milestone 2"
hascode = true
rss = "Description"
rss_title = "Milestone 2"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 2 - Heat Capacity and Time Dependence

\toc

We are currently considering the simple algebraic equilibrium EBM,
$$
(A(CO_2) + B\,T) 4\,\pi\,R_E^2 = S_0\,(1-\alpha)\,\pi\,R_E^2.
$$

As our next step, we want to include a temporal component, i.e., the possibility that temperature $T$ changes in time. Following the general derivations of the [EBM chapter](/milestone2/milestone2_ebm/), we add the temporal term, which is proportional to the temporal derivative $\frac{\partial T}{\partial t} = T_t$. Hence, we add the temporal term with the constant $C$ to our model
$$
(C\,T_t + A(CO_2) + B\,T) 4\,\pi\,R_E^2 = S_0\,(1-\alpha)\,\pi\,R_E^2
$$
to get
$$
C\,T_t + A(CO_2) + B\,T = \frac{S_0}{4}\,(1-\alpha), 
$$
where analysis of the physical units reveal that the Earth's surface area normalized heat capacity $C$ has the units $[J/m^2/K]$.

@@colbox-blue
**Remark 10:** The term heat capacity hints towards a process of "storing the heat". From our observation (and every day experience) we know that depending on the material, the
storage capabilities of heat can be different: For instance when we compare the heat in air and the heat stored in water, or the heat stored in the solid walls of our appartment rooms. For our
full EBM model we want to take this into account and want to have different values of $C$ depending on the location, i.e., if the location is on land, ocean, snow covered
or with sea ice. We want to make it dependent on the geography (from milestone 1).
@@

In general, we need an estimate of the mass density $\rho$ with units $[kg/m^3]$, the specific heat capacity $C_p$ with units $[J/kg/K]$, and the estimated height $d$ in meters $[m]$ of the material column at a given location $x$. 

For the **atmosphere** (air), we have $\rho_{atm} = 1.225$, $C_{p,atm}=1000$, and $d_{atm}=3850$ (here, we take into account that the density of the atmosphere reduces the higher we get and compute an corresponding integrated average height), which gives the value 
$$
\widetilde{C}_{atm} = \rho_{atm}\,C_{p,atm}\,d_{atm} = 4.71625\cdot 10^{6}\,\,[J/m^2/K].
$$

For the **ocean** (mixed layers of water), we have $\rho_{mixed} = 1030$, $C_{p,mixed}=4000$, and $d_{mixed}=70$, which gives the value 
$$
\widetilde{C}_{mixed} = 2.884\cdot 10^{8}\,\,[J/m^2/K].
$$

For the **sea ice**, we have $\rho_{ice} = 917$, $C_{p,ice}=2000$, and $d_{ice}=1.5$, which gives the value 
$$
\widetilde{C}_{ice} = 2.751\cdot 10^{6}\,\,[J/m^2/K].
$$

For the **snow cover**, we have $\rho_{snow} = 400$, $C_{p,snow}=880$, and $d_{snow}=0.5$, which gives the value 
$$
\widetilde{C}_{snow} = 1.76\cdot 10^{5}\,\,[J/m^2/K].
$$

For the **soil**, we have $\rho_{soil} = 1350$, $C_{p,soil}=750$, and $d_{soil}=1$, which gives the value 
$$
\widetilde{C}_{soil} = 1.0125\cdot 10^{6}\,\,[J/m^2/K].
$$

@@colbox-blue
**Remark 11:** These numbers are directly from the paper of Zhuang et al. (2017, Table 1, note that there is a typo in the table for the ocean). However, there is a range of possible numbers available in the literature, depending on the authors.
@@

With these estimates of the heat capacities for different materials, we are almost ready to define the heat capacity values for the geography map from milestone 1. 

But first we introduce a normalization with respect to our time scale that we will need later on for our full EBM. The typical units in physics for time $t$ is seconds $[s]$. However out of convenience we will instead use the unit year $[yr]$ instead. This affects the time derivative of temperature $T_t$ that we need to scale properly with the number of seconds per year. We choose $sec\_per\_year = 3.15576\cdot 10^7$ following the paper by Zhuang et al. (2017), which corresponds to $365,25$ days in the year. To not carry around this scaling factor explicitly, out of convenience, we absorb this parameter into the heat capacity $C(x)$ as it gets multiplied by $T_t$ anyways. 

We hence get, for the different geopgraphy types in the map: 

Land/soil ($geo=1$): $C = (\widetilde{C}_{soil} + \widetilde{C}_{atm})/sec\_per\_year$.

Sea ice ($geo=2$): $C = (\widetilde{C}_{ice} + \widetilde{C}_{atm})/sec\_per\_year$.

Snow cover ($geo=3$): $C = (\widetilde{C}_{snow} + \widetilde{C}_{atm})/sec\_per\_year$.

Lakes/inland sea ($geo=4$): $C = (\widetilde{C}_{mixed}/3 + \widetilde{C}_{atm})/sec\_per\_year$.

Ocean ($geo=5$): $C = (\widetilde{C}_{mixed} + \widetilde{C}_{atm})/sec\_per\_year$.

@@colbox-blue
**Remark 12:** As mentioned, we follow closely the work of Zhuang et al. (2017) and also use the geography data from this work. This geography data does not contain the case $geo = 4$, but defines every water body as ocean ($geo = 5$). It would be interesting to investigate the impact of different bodies of water in the geography. Of course, this only makes sense, when the grid size is small enough to even resolve local (bigger) lakes. 
@@

In summary, we have now our first time dependent (still simplified) EBM
$$
C(x)\,T_t(x,t) + A(CO_2) + B\,T(x,t) = \frac{S_0}{4}\,(1-\alpha).
$$

