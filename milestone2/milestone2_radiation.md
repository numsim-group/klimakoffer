+++
title = "Milestone 2"
hascode = true
rss = "Description"
rss_title = "Milestone 2"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 2 - Radiation modeling

\toc

## Introduction

The energy balance of Earth is strongly impacted by radiation. The following figure shows the energy fluxes in the global Earth atmosphere system.

\fig{/assets/milestone2/ClimateSystem.jpg}
* Source: https://web.archive.org/web/20140421050855/http://science-edu.larc.nasa.gov/energy_budget/ quoting Loeb et al., J. Clim 2009 & Trenberth et al, BAMS 2009

@@colbox-blue
**Remark:** Note that this figure is from 2010. By now the net absorbed energy has risen to $1.0\; W/m^2$. ([Loeb et al. (2021). Satellite and Ocean Data Reveal Marked Increase in Earth’s Heating Rate. Geophysical Research Letters, 48(13)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021GL093047))
@@

In our EBM, we aim to include four effects: 

(i) Incoming solar radiation. 

(ii) Reflection of radiation by the surface and cloud cover. 

(iii) Cooling because of outgoing longwave radiation. 

(iv) Effect of greenhouse gases (mainly $CO_2$).

## Outgoing longwave radiation 

The goal in this section is to define a model/parametrization for the source term $S_{OLW}(T,x,t)$. 

But first, to warm up with the topic, we consider as a toy model an idealized black-body Earth, i.e., we assume that the Earth is a _black body_ that emits and absorbs radiation in the infrared spectrum.

The radiation energy per time that is emitted by a black body can be computed by the Stefan-Boltzmann law of physics 
$$
I = \sigma_{SB}\,T_R^4,
$$
where $I$ is the radiation energy per time per area with units $[W/m^2]$, $T_R$ is the radiation temperature with physical units Kelvin $[K]$, $\sigma_{SB} = 0.56687\cdot 10^{-8}$ is the Stefan-Boltzmann constant with units $[W/m^2/K^4]$.

For this idealized black-body Earth, the temperature is determined when the outgoing radiation is in equilibrium/balance with the incoming stellar radiation. The amount of incoming energy can be roughly estimated as $S_0\,(1-\alpha)\,\pi\,R_E^2$, where $S_0 = 1360 [W/m^2]$ is the solar constant (the mean solar elecromagnetic radiation received on Earth), $\alpha$ is the surface albedo (the amount of the solar radiation that is reflected back to space) with the planetary average being about $\alpha=0.3$ (more details in the [Albedo section](/milestone2/milestone2_albedo/)), and $R_E=6.378\cdot 10^{6} [m]$ is the radius of Earth. 
Note that the solar constant is scaled with the _effective_ area in which the solar radiation is applied on earth: $\pi\,R_E^2$.

Taking into account that the surface of the sphere is $4\pi R_E^2$, the amount of outgoing energy is $\sigma_{SB} T_R^4 4\,\pi\,R_E^2$, and we get our first (simplest version of an) EBM
$$\label{eq:blackbody}
\sigma_{SB} T_R^4 4\,\pi\,R_E^2 = S_0\,(1-\alpha)\,\pi\,R_E^2.
$$
We can relate this simple EBM \eqref{eq:blackbody} to the [general EBM form](/milestone2/milestone2_ebm/#eqebm) by making the assumption of thermodynamic equilibrium (no temporal change $\partial T / \partial t= 0$), no heat diffusion ($d = 0$), and the following choice of source terms: 
\begin{align}
S_{sol} &= S_0\,(1-\alpha)\,\pi\,R_E^2, \\
S_{OLW} &= - \sigma_{SB} T_R^4 4\,\pi\,R_E^2.
\end{align}

We can directly solve \eqref{eq:blackbody} to get the black-body equilibrium temperature
$$
T_R = \left(\frac{S_0}{4}\frac{(1-\alpha)}{\sigma_{SB}}\right)^{\frac{1}{4}}\approx 255\, [K] = -18\, [^\circ C].
$$
This model is indeed as simple as it gets and it is no surprise that the quality of the prediction of an average Earth temperature is quite off. The black-body radiation temperature of Earth would be only about $-18$ degree Celsius, hence, some major modeling improvements are necessary. 

## Budyko's Empirical infrared Model
Budyko (1968) suggested an empirical linear model for the outgoing longwave radiation 
$$
I \sim A + B\, T.
$$
It is in general motivated by available observational data, shown in the next figure

\fig{/assets/milestone2/OutgoingLongwaveRadiation.png}
* Figure is from the book _North, G. R., & Kim, K. Y. (2017). Energy balance climate models. John Wiley & Sons_ and is originally from [Graves, C. E., Lee, W. H., & North, G. R. (1993). New parameterizations and sensitivities for simple climate models. Journal of Geophysical Research: Atmospheres, 98(D3), 5025-5036](https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/92JD02666?casa_token=X0WG_pxk8AUAAAAA:2mvPv6HgmsA467qq44RYKY8WrJZLh_Bl-lN2kzgdBLJi3-xSVh0il6g-p1PSlxda51H8YVdkx1dsxSI).

The figure shows infrared radiation density plots averaged monthly, measured by satellite compared to the surface temperature at the same month and location. 
(a) shows the whole sky (including clouds) and (b) shows only the clear (cloudless) sky.

If we consider the temperature in units Kelvin, we can fit the observed data with the linear model by choosing good constants $A$ and $B$ to get 
$$
I_{IR/OLW} = A + B\,(T - 273.15),
$$
with $A=210.3$ as the radiative cooling in units $[W/m^2]$, and $B=2.15$ the radiative cooling feedback with units $[W/m^2/K]$. 
It is important to note that the choice of these parameters has a direct impact on the outgoing radiation and hence on the cooling. Several others have fitted the data differently, hence some range of choices for $A$ and $B$ is available. The values we select are from the paper by Zhuang et al. (2017).

We are now able to consider a second, but hopefully improved toy EBM. We replace the crude black-body radiation with a phenomenological approximation of the outgoing radiation (the Budyko model) to get 
$$
(A + B\,(T_R - 273.15)) 4\,\pi\,R_E^2 = S_0\,(1-\alpha)\,\pi\,R_E^2,
$$
which can be directly solved again to get the equilibrium temperature
$$
T_R = \frac{\frac{S_0}{4}(1 - \alpha) - A}{B} + 273.15 \approx 13\,[^\circ C].
$$

@@colbox-blue
**Remark 8:** As discussed, we have introduced a data driven model for our outgoing longwave radiation energy change. We introduced an ansatz with parameters (parametrization) and used real world measurements/observations to select the parameters, i.e., to fit the model.
@@

@@colbox-blue
**Remark 9:** The phenomenological model of Budyko drastically improves the temperature prediction of Earth from the crude black-body temperature of $-18$ to the temperature $13^\circ C$, compared to the observed value of about $14.5^\circ C$.
@@

## Effect of greenhouse gases ($CO_2$)

For the interested reader, we refer to chapter 4 of the book by Kim and North, "Energy Balance Climate Model", (2017, Wiley) and the following paper 

> [Myhre, G., Highwood, E. J., Shine, K. P., & Stordal, F. (1998). New estimates of radiative forcing due to well mixed greenhouse gases. Geophysical research letters, 25(14), 2715-2718](https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/98GL01908).

The Earth's atmosphere is composed of several gases, with the most abundant ones being nitrogen ($N_2$), oxygen ($O_2$), and argon ($Ar$). These gases do not strongly absorb infrared radiation, which is important because the Earth's surface emits infrared radiation as it cools down after being heated by the sun.
$N_2$ and $O_2$ are diatomic molecules, meaning they consist of two atoms chemically bonded together. Diatomic molecules have no permanent dipole, which means they have no separation of electric charge and therefore do not strongly interact with infrared radiation.
Argon, on the other hand, is a monatomic gas, meaning it consists of individual atoms rather than molecules. However, it does not absorb infrared radiation because it does not have any modes of rotation or vibration in the infrared spectrum.

Overall, the lack of strong absorption of infrared radiation by these main constituents of the atmosphere allows heat to escape from the Earth's surface and be radiated out to space, helping to regulate the planet's temperature.

$H_2O$ molecules, on the other hand, have a permanent dipole moment because the distribution of electrons in the molecule is asymmetric. As a result, $H_2O$ molecules respond strongly to passing electromagnetic waves, including those in the infrared spectrum, which leads to the absorption of infrared radiation. 

Molecules such as $CO_2$, $CH_4$, and $NO_2$ do not have a permanent dipole moment because they are symmetrical in shape. However, they can still absorb infrared radiation through the phenomenon of induced dipole moments. When an infrared photon passes near one of these molecules, it can cause the electrons in the molecule to shift slightly, resulting in a temporary dipole moment. This temporary dipole moment can then interact with the passing electromagnetic wave, leading to the absorption of infrared radiation.

The ability of these molecules to absorb infrared radiation is significant because it allows them to contribute to the greenhouse effect. The greenhouse effect is the process by which certain gases in the Earth's atmosphere, including $CO_2$, $CH_4$, and $NO_2$, trap heat and warm the planet's surface. Without this natural process, the Earth's average temperature would be much lower and life as we know it would not be possible. However, when these gases are present in excess, they can cause an imbalance in the greenhouse effect, leading to climate change.

@@colbox-blue
**Comment:** The effect of greenhouse gases is very complex and there is a lot of ongoing research. $H_2O$ in particular is very important but ever so complex due to it being dependent on many other components and processes, i.e., the full water cycle. We will consider in our EBM only the effect of $CO_2$ on the radiation. Infrared absorbtion means that the efficiency of our outgoing radiation (cooling effect) decreases (which causes a higher equilibrium temperature). This decrease in efficiency is topic of many investigations and studies, e.g., by the IPCC.
@@

In this course, we follow the paper by Myhre et al. (1998) to define our parametrization of the greenhouse gas effect. From this paper, we first look at the effect of the amount of greenhouse gases, as plotted in the following figures:

\fig{/assets/milestone2/CO2_forcing.png}
\fig{/assets/milestone2/CH4_forcing.png}
\fig{/assets/milestone2/N2O_forcing.png}
 * Figures were generated with data from Myhre et al. (1998).
 
The figures show radiative forcing as a function
of concentration in $[ppmv]$ ("parts per million volume") for $CO_2$ and in $[ppbv]$ ("parts per billion volume") for $CH_4$ and $N_2O$.

The following table gathers simplified expressions to compute (fit) the radiative forcing caused by different greenhouse gases using data from Myhre et al. (1998) and the function
\begin{align}\label{eq:ffunciton}
f(M,N) = 0.47 \ln \left( 1 + 2.01 \times 10^{-5} (M \, N)^{3/4} + 5.31 \times 10^{-15} M (M \, N)^{1.52} \right).
\end{align}

* **Table:** Simplified expressions for the radiative forcing ($\Delta F$) in $W/m^2$ with coefficients of the IPCC report (1990) and Myhre et al. (1998). In the expressions, $C$ is CO$_2$ concentration in $[ppmv]$, $M$ is CH$_4$ concentration in $[ppbv]$, $N$ is N$_2$O concentration in $[ppbv]$, and $X$ is Chlorofluorocarbons (CFCs) concentration in $[ppbv]$. The subscript $0$ denotes unperturbed (reference) concentrations. The function $f$ is given in \eqref{eq:ffunciton}. Adapted from Myhre et al. (1998).
\begin{align*}
\begin{array}{cccc}
\hline
\text{Trace gas}   & \text{Radiative forcing $\Delta F$}     & \alpha_{\text{IPCC}} & \alpha_{\text{Myhre}} \\
\hline
\text{CO}_2  &   \alpha \ln (C/C_0) & 6.3 & 5.35      \\
\text{CH}_4  &   \alpha \left(\sqrt{M}-\sqrt{M_0} \right) - \left(f\left(M,N_0\right)- f\left(M_0,N_0\right) \right) & 0.0036        & 0.0036   \\
\text{N}_2\text{O} &  \alpha \left(\sqrt{N}-\sqrt{N_0} \right) - \left(f\left(M_0,N\right)- f\left(M_0 ,N_0\right) \right)  &0.14 & 0.12 \\
\text{CFC-11} & \alpha(X-X_0) & 0.22 & 0.25\\ 
\text{CFC-12} & \alpha(X-X_0) & 0.28 & 0.33\\ 
\hline
\end{array}
\end{align*}



As mentioned, we only consider the effect of $CO_2$ in our model and hence choose the approximation 
$$
\Delta F = \alpha_{\text{Myhre}}\, \ln(CO_2/CO_2(t_0)),
$$
where 
$$
\alpha_{\text{Myhre}} = 5.35\,[W/m^2],
$$
$CO_2$ is the concentration in $[ppm]$, and $CO_2(t_0) = 315\, [ppm]$ is the reference concentration in the year $t_0 = 1950$.

In summary, we get the following parametrization of our outgoing longwave radiation source term 
$$
S_{OLW}(T,x,t,CO_2) = -(A + B\, (T - 273.15) - \Delta F(CO_2)),
$$
which we reformulate into the shorthand notation
$$
S_{OLW}(T,x,t,CO_2) = - (A(CO_2) + B\,T),
$$
where
$$
A(CO_2) := 210.3 - 5.35\,\ln(CO_2/315))\,\,[W/m^2],
$$
and 
$$
B = 2.15\,\,[W/m^2/K],
$$
where we further made the assumption that the unit of the temperature $T$ is  $[^\circ C]$ instead of Kelvin.

