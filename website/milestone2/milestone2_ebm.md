+++
title = "Milestone 2"
hascode = true
rss = "Description"
rss_title = "Milestone 2"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 2 - The Energy Balance Model

\toc

## Introduction

The Energy Balance Climate Models (or Energy Balance Models (EBM) in short) were
introduced in the 1960s by Budyko and Sellers independently (with some variations). They are thus sometimes referred to as Budyko-Sellers equations (or model). 

EBMs where mainly used up to the 1970s, but later replaced more and more by the GCMs. Nowadays EBMs belong to the class of simple climate models, but
are still in use by some researchers to date, with publications in 2020s. Interestingly enough, EBM type models are used to estimate climates of (habitable) exoplanets (where naturally not many detailed data is known). Their disadvantage (being too simple) is sometimes an advantage when trying to analyse the effect of single processes or at least
get an intuition about how they act. To quote the book by Kim and North (2017) "Energy balance climate models", Wiley:

> In some cases such as perturbations of the surface temperature field due to small changes in greenhouse gases or the Earth's orbital elements, they can be surprisingly helpful even to a quantitative extend.

We will learn that EBMs (and to some extend all other models as well) depend on empirical parameters. The strength
of the EBM is that the model's use of phenomenological coefficients allows it to be grounded on and close to large scale observations. The strong disadvantage
of this approach is that by adjusting parameters to fit current observations, we do not really know how these coefficients change when the climate changes. We
have somewhat **departured from first principles**.

@@colbox-blue
**Remark 6:** It is interesting to realise that "departure from first principles" is currently on the rise in the context of purely data-driven modeling based
on methods from machine-learning and artificial intelligence. There is currently a strong research push to enhance/generate parametrizations and even full sub-models/components
of GCMs/ESMs that are data-driven. As it is well known, adding more and more components and physical processes also adds more parameters (and with this improved realism), 
but also the tendency (or at least danger) of over fitting to existing observational data, which would diminish the predictive power of the models for scenarios they are not trained in.
@@

## Energy Balance Law

As we have discussed before, a balance law has the general form 
$$
\frac{\partial u}{\partial t} + \vec{\nabla}\cdot\vec{f} = S,
$$
in a spatial domain $\Omega$, with boundary conditions for $u$ at the surface $\partial\Omega$ and initial conditions $u(x,t=0) = u_0(x)$.

The modeling process involves the choice of $u$, $\vec{f}$, and $S$. As the name EBM suggests, we consider as our unknown the
energy, more precisely the internal energy. The internal energy of a body/fluid is proportional to its temperature. In the context of climate modelling we are interested
in simulating said temperature, more precisely we are interested to approximate the **surface temperature** $T(x,t)$.

We hence make our first modeling step and choose $u(x,t) \sim T(x,t)$, i.e. 
$$
u(x,t) = C(x)\,T(x,t),
$$
where the parameter $C(x)$ is a scaled **heat capacity**, which combined with temperature gives the model of internal energy. 

@@colbox-blue
**Important note:** We use the variable $x$ as a placeholder notation for the coordinate variables. We want to solve the EBM on the sphere surface. Hence, $x$ refers to the spherical coordinates (latitude/colatitude and longitude).
@@

As mentioned, we are interested in the (surface) temperature $T(x,t)$ and its temporal evolution - the fluxes can be modeled analogously to the fluxes of the heat transfer equation $\vec{f}\sim-\vec{\nabla} T$, which leads to the model of the EBM flux as 
$$
\vec{f} = -D(x)\,\vec{\nabla} T,
$$
where the diffusion coefficient has positive entries, $D(x)\in\mathbb{R}_+$.

In this course, we consider an advanced version of the EBM where detailed solar/stellar radiation is the incoming energy source term and the outgoing longwave radiation (in the infrared) is an energy sink term, i.e., the source term has two major parts
$$
S(u,x,t) = S_{OLW}(T,x,t) + S_{sol}(x,t).
$$

With these first modeling steps/decisions, the general form of the EBM is 
$$\label{eq:EBM}
C(x)\,\frac{\partial T}{\partial t} - \vec{\nabla}\cdot(D(x)\,\vec{\nabla} T) = S_{OLW}(T,x,t) + S_{sol}(x,t).
$$

@@colbox-blue
**Remark 7:** We already noted that we consider the EBM on the sphere, i.e., on the Earth's surface. This means that for the
discretization we will need to get more into the detail on how the differential operators $\vec{\nabla}$ are defined in spherical coordinates. 
@@

The modeling process is not finished yet, as we need detailed definitions of the heat capacity $C(x)$, the diffusion coefficient $D(x)$, and the sources $S_{OLW}(T,x,t)$, $S_{sol}(x,t)$. 

We will focus next on all of these sub-models and parametrizations, except for the diffusion coefficient $D(x)$ and the approximation of the diffusion operator. This will be discussed in [milestone 5](/milestone5/milestone5_menu/) in detail.
