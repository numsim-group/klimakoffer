+++
title = "Milestone 5"
hascode = true
rss = "Description"
rss_title = "Milestone 5"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 5 - Planetary Heat Transfer

\toc

## Transport in the Atmosphere and Oceans

As a next improvement of our EBM, we want to connect the temperatures of the each individual grid cell to form a globally coupled system that models a smooth surface temperature field of the whole Earth. If we consider a fully coupled system where all grid cells are connected (and not independently solved as tested in [milestone 4](/milestone4/milestone4_menu/)), heat will get transported from warm regions to colder regions, i.e., heat will be transported polewards. 
The [distribution of the solar forcing term](/milestone2/milestone2_results/#solar_forcing_animation) clearly shows that there is more net solar radiation in tropical latitudes than in the polar regions. Thus, some heat transfer polewards is expected. On Earth, the atmosphere and the oceans include mechanisms and processes that cause such a poleward heat transport. For the interested reader, we refer for instance to the following introductory videos that discuss some of the important processes.
 
* Global atmospheric circulation on Earth
~~~
<iframe width="560" height="315" src="https://www.youtube.com/embed/xqM83_og1Fc" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
~~~
* Effect of Coreolis force in large-scale vortices
~~~
<iframe width="560" height="315" src="https://www.youtube.com/embed/PDEcAxfSYaI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
<br>
<br>
~~~

An important property of moving fluid flows is turbulent motion. Most flows in nature (and engineering) are turbulent, which means multiscale in space and time. From very small eddies in the range of milimeters up to large vortical structures of size hundreds or thousands of meters. Turbulent flows are highly sensitive to small scale disturbances (this is known as the butterfly effect), almost chaotic and thus hard to predict. One of the very important mechanisms of turbulent flows is increased mixing efficiency of cold and warm fluids, i.e., an enhanced heat transfer and heat distribution. This formally requires the numerical resolution (very small grid cells) of all small and large eddies. Small grid cells of course make the simulations very computationally intense and unfortunately for realistic flow scenarios too expensive, even on today's largest super computers. Hence, the full resolution of turbulent fluids is an unsolved problem in fluid mechanics. The current strategy in computational fluid dynamics (and also atmospheric science) to deal with this problem is to introduce additional parametrizations (also often coined turbulence models) that account for the effects of small eddies that cannot be resolved by the grid, so-called subgrid effects. Hence, all modern ESM include a full atmosphere simulation model, which itself includes parametrization for subgrid scale turbulence. 

@@colbox-red
**Bad news:** Unfortunately, a full atmospheric model is out of the scope of the course.
@@

@@colbox-green
**Good news:** It is possible to model the increased mixing of the turbulent fluid and the general convection transport of heat towards the pole by a simplified (but highly parametrized) process, namely by a diffusion process.
@@

## Modeling Turbulent and Polewards Heat Transport with Diffusion

In his pioneering work, Budyko 
> [Budyko, M. I. (1969). The effect of solar radiation variations on the climate of the Earth. tellus, 21(5), 611-619.](https://www.tandfonline.com/doi/pdf/10.3402/tellusa.v21i5.10109)
considered a simple 1D (latitude dependent) EBM that included a simple term that models heat transfer accross latitude zones. From satellite data one can reconstruct that the main heat transfer depends on the latitude and is polewards. The next figure, for instance, shows the heat transfer in the northern hemisphere. It is also worth noting that the heat transfer in the atmosphere is different compared to the ocean. 

\fig{/assets/milestone5/heat_transfer_north.png}
* Heat transfer in the northern hemisphere. Generated with data from "Maslin (2013). 'Climate: A Very Short Introduction'"

The idea of Budyko was to include an algebraic term that is proportional to 
$$
\sim D\,[T(\lat) - T_{avg}],
$$
where $T(\lat)$ is the temperature at latitude $\lat$, $T_{avg}$ is the mean global surface temperature, and $D$ is a parameter that can be adjusted. The next figure shows a sketch of the model with energy fluxes across the latitude zones (rings), where the coefficient $D$ is denoted with $C$ in this figure.  

\fig{/assets/milestone5/heat_transfer_budyko_sellers.png}
* Figure from [https://pages.jh.edu/lhinnov1/paleoguide/tutorial2.html](https://pages.jh.edu/lhinnov1/paleoguide/tutorial2.html).


Sellers,
>[Sellers, W. D. (1969). A global climatic model based on the energy balance of the earth-atmosphere system. Journal of Applied Meteorology and Climatology, 8(3), 392-400.](https://journals.ametsoc.org/downloadpdf/journals/apme/8/3/1520-0450_1969_008_0392_agcmbo_2_0_co_2.pdf)
proposed a model for the heat transport that is based on a diffusion operator. Diffusion or heat conduction is a model where the heat flux is proportional to the gradient of the temperature 
$$
\sim D\,\Nabla T,
$$
and the diffusion operator is of the form 
$$
\sim \Nabla\cdot (D\,\Nabla T),
$$
where $D$ is now the so-called diffusion coefficient. The diffusion coefficient needs to include the effect of oceanic heat transport, turbulent mixing, the sensible heat transport and the latent heat transport. Sensible heat transfer occurs when objects with different temperatures are in direct contact (e.g. surface of the land with the atmosphere) and exchange heat. Latent heat cannot be measured directly, but is an important form of heat transfer that is relevant when a material changes its state (from liquid to gas or from ice to liquid). 

@@colbox-blue
**Remark:**
It is clear that all these complex processes cannot be phenomenologically modeled by a single parameter $D$. Coming back to our discussion on the different modeling strategies, we already had models based on first principles. These models are (supposed to be) universally valid and based on fundamental laws of nature with no tunable parameters. A next step was to approximate processes that can be measured with simplified models. These empirical models are tuned with parameters to approximate the behavior of the processes observed in measurements. The empircal models are still grounded in reality as they are calibrated with measurements, but include some for of parameter tuning. As an example, we had the radiation model of Budyko in milestone 2.

Due to the complexity of the heat transport and the many mechanisms and processes involved, we do not consider first principle modeling, nor empirical models for the heat transfer in our EBM. 
Instead, we resort now to another strategy in modeling. We add a simple replacement model, which does not directly mimic the processes, but mimics the effect of the processes on the solution of the EBM. To be precise, we aim to model the effect of the poleward heat transfer on the resulting EBM temperature. This simple replacement model has again parameters, which cannot be tuned to measurements and observations of processes. Instead, the idea is to tune the parameters of the model such, that the resulting solution of the model (the surface temperature) behaves as realistic as possible.  

It is important to understand that this type of modeling is driven by data of the desired solution (temperature distribution), i.e., we fit the model to existing solutions (for instance from observations) as good as possible. So we do not tune the model of the process, but tune the outcome of the whole model. It is thus clear that such type of modeling is highly heuristic and should be used with extreme caution, as its predictive power is hard to gauge, i.e., expert knowledge and validation is necessary to give such an approach credibility. 
@@


## Diffusion Operator in Spherical Coordinates

As we solve our EBM on the surface of Earth, i.e., in spherical coordinates, we need to define the gradient operators and the divergence operators in our natural coordinate system. 

The gradient and divergence operators defined in Cartesian coordiates read as,
\begin{align}\label{eq:opsCartesian}
\Nabla T &:= \partialderiv{T}{x} \hat{x} + \partialderiv{T}{y} \hat{y} +  \partialderiv{T}{z} \hat{z},
\\
\Nabla \cdot \vec{F} &:= \partialderiv{F_x}{x} + \partialderiv{F_y}{y} + \partialderiv{F_z}{z},
\end{align}
where $T$ is a scalar field, $\hat{x}$, $\hat{y}$, and $\hat{z}$ are the unit vectors in the $x$, $y$, and $z$ directions, respectively, and $\vec{F}=F_x \hat{x} + F_y \hat{y} + F_z \hat{z} = (F_x, F_y, F_z)$ is a vector field.

We can transform the gradient and divergence operators to the spherical coordinate system by computing the relation between the unit vectors and partial derivatives of both coordinate systems.

For instance, for the unit vectors we have
<!-- \begin{align*}
\hat{r} &= \sin \colat \cos \long \hat{x} + \sin \colat \sin \long \hat{y} + \cos \colat \hat{z}
\\
\hat{\colat} &= \cos \colat \cos \long \hat x + \cos \colat \sin \long \hat{y} - \sin \colat \hat z
\\
\hat{\long} &= -\sin \long \hat x+ \cos \long \hat y
\end{align*} -->
\begin{align}\label{eq:unitvectors}
\hat{x} &= \sin \colat \cos \long \hat{r} + \cos \colat \cos \long \hat{\colat} - \sin \long \hat{\long}
\\
\hat{y} &= \sin \colat \sin \long \hat{r} + \cos \colat \sin \long \hat{\colat} + \cos \long \hat{\long}
\\
\hat{z} &= \cos \colat \hat{r}- \sin \colat \hat{\colat},
\end{align}
and for the partial derivatives we have
\begin{align}\label{eq:partialderivs}
\partialderiv{}{x} &= \sin \colat \cos \long \partialderiv{}{r} + \frac{\cos \colat \cos \long}{r} \partialderiv{}{\colat} - \frac{\sin \long}{r \sin \colat} \partialderiv{}{\long}
\\
\partialderiv{}{y} &= \sin \colat \sin \long \partialderiv{}{r} + \frac{\cos \colat \sin \long}{r} \partialderiv{}{\colat} + \frac{\cos \long}{r \sin \colat} \partialderiv{}{\long}
\\
\partialderiv{}{z} &= \cos \colat \partialderiv{}{r} -  \frac{\sin \colat}{r} \partialderiv{}{\colat}.
\end{align}

Inserting \eqref{eq:unitvectors} and \eqref{eq:partialderivs} into \eqref{eq:opsCartesian}, we obtain the gradient and divergence operators in spherical coordinates,
\begin{align}
\Nabla T &= \frac{\partial T}{\partial r} \hat{r} + \frac{1}{r\sin(\colat)}\frac{\partial T}{\partial \long} \hat{\long} + \frac{1}{r}\frac{\partial T}{\partial \colat} \hat{\colat}
\\
\Nabla \cdot \vec{F} &= 
\frac{1}{r^{2}}\partialderiv{}{r}({F}_r r^{2})
+
\frac{1}{r\sin(\colat)} \frac{\partial {F}_{\long}}{\partial \long}
+ 
\frac{1}{r\sin(\colat)}\frac{\partial}{\partial \colat}  ({F}_{\colat}\sin(\colat))
\end{align}


The gradient operator is also commonly written in the vector form $(r,\long,\colat)$ as
\begin{align}
\label{eq:12}
\Nabla T = \left(\frac{\partial T}{\partial r}, \frac{1}{r\sin(\colat)}\frac{\partial T}{\partial \long}, \frac{1}{r}\frac{\partial T}{\partial \colat}\right),
\end{align}


We obtain the diffusion operator in spherical coordinates by combining the expressions for the gradient and divergence operators to obtain
\begin{align}
\Nabla \cdot (D\Nabla T) = 
\frac{1}{r^2} \partialderiv{}{r} \left( D r^2 \partialderiv{T}{r} \right)
+
\frac{\csc^{2}(\colat)}{r^2} \frac{\partial}{\partial \long}\Bigl(D\frac{\partial T}{\partial \long}\Bigr) 
+ 
\frac{\csc (\colat)}{r^2}
\frac{\partial}{\partial \colat}\left(D \sin (\colat)\frac{\partial T}{\partial \colat}\right)
\end{align}

Note that these derivations are in full three-dimensional space. However, in our EBM, we consider only the surface of the Earth, i.e. a $2D$ approximation of the surface temperature with the choice of the radius coordinate $r = R_E$. Hence, heat transfer in vertical direction (along the radial coordinate $r$) is neglected and only the variations in latitude and longitude directions are considered. Furthermore, we absorb the scaling with the Earth radius $R_E$ into the diffusion coefficient and simplify the expression to obtain
\begin{align}\label{eq:diffterm}
\Nabla \cdot (D\Nabla T) = 
\underbrace{
    \csc^{2}(\colat) \frac{\partial}{\partial \long}\Bigl(\diffcoeff\frac{\partial T}{\partial \long}\Bigr) 
}_{\text{Term}~1}
+ 
\underbrace{
    \frac{\partial}{\partial \colat}\Bigl(\diffcoeff\frac{\partial T}{\partial \colat}\Bigr)
}_{\text{Term}~2}
+ 
\underbrace{
    \cot(\colat)\diffcoeff\frac{\partial T}{\partial \colat}.
}_{\text{Term}~3}
\end{align}

@@colbox-blue
**Remark:** The diffusion coefficient that we use in our $2D$ EBM model is scaled with the Earth radius: $\diffcoeff := D / R_E^2$.
@@


## Choice of the Diffusion Coefficient in the 2D EBM

Our final form of the $2D$ EBM with heat diffusion reads as 
$$
C(x) \partialderiv{T}{t} + A(CO_2) + B T - \Nabla \cdot (D\Nabla T) = S_{sol}(x,t).
$$
respectively, in spherical coordinates
\begin{align}\label{eq:ebm_spherical}
C(\colat,\long) \partialderiv{T}{t} &+ A(CO_2) + B T 
\\
 &- \left[
    \csc^{2}(\colat) \frac{\partial}{\partial \long}\Bigl(\diffcoeff\frac{\partial T}{\partial \long}\Bigr) 
+ 
    \frac{\partial}{\partial \colat}\Bigl(\diffcoeff\frac{\partial T}{\partial \colat}\Bigr)
+ 
    \cot(\colat)\diffcoeff\frac{\partial T}{\partial \colat}\right]
= S_{sol}(\colat,\long,t).
\end{align}

@@colbox-blue
**Remark:** An analysis of the physical dimension of the diffusion coefficient shows that $D$ has the SI units $\left[\frac{W}{K}\right]$ (Watts per Kelvin), and $\diffcoeff$ has the SI units $\left[\frac{W}{m^2 K}\right]$ (Watts per Kelvin per square meter).
@@

As discussed above, the heat diffusion term needs to model the complex mechanism of poleward transport of heat. There are many different choices of $D$ available in literature, from simple global constant values to more complex approximations. Here, again, we follow the paper by Zhuang et al. 
> [Zhuang, K., North, G. R., & Stevens, M. J. (2017). A NetCDF version of the two-dimensional energy balance model based on the full multigrid algorithm. SoftwareX, 6, 198-202.](https://www.sciencedirect.com/science/article/pii/S2352711017300262)

We distinguish between oceanic heat transport and heat transport over land and snow/ice covered areas. We also account for a difference of the heat transfer in the northern and southern hemisphere due to the asymmetric distribution and sizes of the land masses. The following figure shows a sketch of the energy fluxes due to heat transfer accounted in our model (here, the diffusion coefficients are denoted with $K$ instead of $\diffcoeff$) 

\fig{/assets/milestone5/heat_transfer_sketch_ebm.png}
* Figure from [TransEBM v. 1.0: description, tuning, and validation of a transient model of the Earth's energy balance in two dimensions](https://gmd.copernicus.org/articles/14/2843/2021/) by [Elisa Ziegler](https://orcid.org/0000-0002-7252-3332) and [Kira Rehfeld](https://orcid.org/0000-0002-9442-5362), [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Thus, for the actual values of the diffusion coefficients we differentiate between grid cells that represent ocean, and all other types. Furthermore, for non oceanic grid cells, we differentiate if we are in the northern or southern hemisphere. 

The values for oceanic grid cells in physical unit $[W/m^2/K]$ are 
$$
\diffcoeff = \diffcoeff_{\text{ocean,poles}} + (\diffcoeff_{\text{ocean,equ}} - \diffcoeff_{\text{ocean,poles}})\,sin^5(\colat),
$$
with $\diffcoeff_{\text{ocean,poles}} = 0.4$ and $\diffcoeff_{\text{ocean,equ}} = 0.65$. 

The values for non oceanic grid cells in the northern hemisphere in physical units $[W/m^2/K]$ are
$$
\diffcoeff = \diffcoeff_{\text{NP}} + (\diffcoeff_{\text{equ}} - \diffcoeff_{\text{NP}})\,sin^5(\colat),
$$
with $\diffcoeff_{\text{NP}} = 0.28$ and $\diffcoeff_{\text{equ}} = 0.65$.

The values for non oceanic grid cells in the southern hemisphere in physical units $[W/m^2/K]$ are
$$
\diffcoeff = \diffcoeff_{\text{SP}} + (\diffcoeff_{\text{equ}} - \diffcoeff_{\text{SP}})\,sin^5(\colat).
$$
with $\diffcoeff_{\text{SP}} = 0.20$ and $\diffcoeff_{\text{equ}} = 0.65$.

@@colbox-blue
**Remark:** We stress again, that the choice of the diffusion coefficients is kind of arbitrary and highly tuned to fit the resulting temperature fields as good as possible. Thus, the choice of a ramp up function between the values at the equator and the values at the poles proportional to $sin^5(\colat)$ has no deeper meaning but is a modeling/tuning choice. The resulting distribution of the diffusion coefficients over ocean and non oceanic grid cells is shown in the following figure 

\fig{/assets/milestone5/diffusion_coefficient.png}

Note that this plot is showing the latitude rather than the colatitude.
@@

