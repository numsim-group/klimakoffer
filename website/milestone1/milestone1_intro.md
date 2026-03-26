+++
title = "Milestone 1"
hascode = false
rss = "Description"
rss_title = "Milestone 1"
rss_pubdate = Date(2019, 5, 1)

tags = ["climatesystem"]
+++

# Milestone 1 - Introduction to Climate
\toc

## Components of the Climate System

The Earth's climate is a complex system that redistributes the energy from the sun throughout the planet. This energy transfer is carried out by various components, including the atmosphere, hydrosphere, cryosphere, land surface, and biosphere. Each of these components plays a crucial role in regulating the Earth's climate, and changes to any one of them can have far-reaching effects on the system as a whole. 
In addition, human activities affect climate on Earth significantly and are nowadays considered an additional component: the anthroposphere.
Understanding how these components interact with each other is essential to comprehending the complex processes that govern the Earth's climate.

\fig{/assets/milestone1/ClimateSystem.png}
<!-- ![](/assets/milestone1/ClimateSystem.png) -->
* Figure from lecture notes: [Stocker, "Introduction to Climate Modeling". Universität Bern](https://climatehomes.unibe.ch/~stocker/papers/stocker18icm.pdf), published by [Springer](https://link.springer.com/book/10.1007/978-3-642-00773-6).

### Atmosphere

The atmosphere is the gaseous layer above the Earth's surface and is composed of various substances in gaseous, liquid (such as water and aerosols), or solid (like dust) forms. It plays a crucial role in regulating the Earth's climate and supporting life on the planet.

Some of the important atmospheric processes include:

* Weather patterns, which involve the short-term changes in temperature, humidity, wind, and precipitation in a particular region.
* Radiation balance, which involves the absorption and reflection of incoming solar radiation, the Earth's albedo, and the emission of radiation from the Earth's surface.
* Cloud and precipitation formation, which are critical processes that help regulate the planet's climate.
* Turbulent atmospheric flow, which plays a crucial role in transferring heat and mixing various atmospheric constituents, thereby influencing weather patterns, cloud formation, and precipitation.

Overall, these processes work together to create a complex and dynamic atmosphere that is vital for sustaining life on Earth.

### Hydrosphere 

The hydrosphere comprises all forms of water on and below the Earth's surface, including oceans, rivers, lakes, and groundwater reservoirs. It also includes the entire global water cycle, which begins when precipitation reaches the surface.

Several important processes occur within the hydrosphere, some of which are:

* Ocean water transport, which is responsible for the movement of water masses and heat around the globe.
* Changes in the inflow of water into different ocean basins, which can have a significant impact on ocean currents and climate patterns.
* The connection between the oceans and the atmosphere, which involves the exchange of water vapor and other gases.
* The absorption of carbon dioxide (CO2) by the oceans, which is the most crucial reservoir for carbon. The oceans absorb more CO2 than the atmosphere and the terrestrial biosphere (plants and animals) combined.


### Cryosphere

The cryosphere refers to all forms of ice in the Earth's climate system, including ice masses, ice shelves, sea ice, glaciers, and permafrost. The amount of ice present in the cryosphere has a significant impact on the hydrosphere, as it serves as a long-term water reserve that can affect water availability in different regions of the world. In addition, the cryosphere plays a crucial role in regulating the Earth's radiation balance through its effect on the planet's albedo, or the reflection of incoming solar radiation. Changes in the cryosphere can lead to alterations in the planet's albedo, which can have significant consequences for global temperature and climate patterns. Therefore, understanding the dynamics of the cryosphere and its relationship with other components of the climate system is essential for predicting and mitigating the effects of climate change.


### Land Surface 

The Land Surface refers to the solid portion of the Earth's crust. The location and positioning of the continents have a significant impact on the distribution of climatic zones across the planet, as well as the formation and direction of ocean currents. Additionally, the Land Surface plays a critical role in regulating the Earth's radiation balance through its effect on the planet's albedo, which can vary depending on the type of terrain (e.g., sand versus rocks). The Land Surface also acts as a reservoir of dust particles that can interact with the atmosphere and affect weather patterns. Understanding the Land Surface and its interactions with the other components of the Earth's climate system is essential for predicting and mitigating the effects of climate change.

### Biosphere 
The Biosphere refers to the organic cover of the Earth's land masses, including vegetation, soil, and marine organisms. It has a strong impact on carbon exchange between different parts of the Earth and can significantly affect the concentration of CO2 in the atmosphere. Some of the key processes that the Biosphere is involved in include:

* Changing the reflectivity (Albedo) of the Earth's surface, which can have a significant impact on the planet's radiation balance.
* Regulating the transfer of water vapor between the land and the atmosphere, which can affect weather patterns and climate.
* Vegetation acts like a rough surface, which can impact atmospheric flow and the exchange of momentum between the land and the atmosphere.

Understanding the role of the Biosphere in the Earth's climate system is critical for predicting and mitigating the effects of climate change.

### Anthroposphere

The Anthroposphere is the term used to describe all the interactions and activities of humans that change existing processes or create new ones within the Earth's climate system. This includes activities such as the high-rate emission of substances, changes in land use, including deforestation, desertification, the conversion of natural habitats into constructed areas, and the drainage of marshes. These human activities can have significant impacts on the Earth's climate and ecosystems, affecting factors such as greenhouse gas emissions, water availability, and the balance of biodiversity. 

Apart from the components mentioned earlier, there are additional components that can significantly impact the Earth's climate system. These components are based on rare and extreme events, such as volcanic eruptions, which can release large amounts of gases and particles into the atmosphere. These emissions can have short-term effects on the planet's temperature and weather patterns by altering the amount of solar radiation that reaches the Earth's surface. Volcanic eruptions can also have long-term impacts by contributing to changes in the Earth's albedo and carbon cycle. Therefore, while such events may be infrequent, they play an important role in shaping the Earth's climate system.

A complete climate model contains all
of the above components. Furthermore,
the interaction between all these components is
necessary (for instance, water vapour exchange).


@@colbox-red
**Bad news**

Every single process requires expert level
of research to understand the physics
and to generate model abstractions. We
could easily fill full lecture courses
on the individual components: the physics, 
the mathematical model, the numerical
algorithms and their implementation, etc.


As a result, a fully coupled global climate model
is out of the scope for this course!
@@

@@colbox-green
**Good news**

Not all research questions in climate science
require the full model. It is, however, part of the
scientific work to select a valid sub-set of
components and processes to get valid
answers that are scientifically robust.
@@

**Comment**: For the interested student, we
recommend to attend other specialized courses
available at UoC or read detailed
lecture notes such as for instance: [Stocker, "Introduction to Climate Modeling". Universität Bern](https://climatehomes.unibe.ch/~stocker/papers/stocker18icm.pdf), published by [Springer (closed access, not available on the university network)](https://link.springer.com/book/10.1007/978-3-642-00773-6).

With this introductory discussion on the climate
system, we come to the conclusion that to achieve
our goal "implement your own climate model from scratch" we need simplifications. To be more
precise: we will need strong simplifications!!


## Hierarchy of Climate Models

The combination of Atmosphere and Hydrosphere models
(in particular oceans) with high-fidelity is termed
Global Climate Model / General Circulation Model (GCM).
As discussed above, these are, however, only two
components out of many. Hence, nowadays
the fully coupled models that includes more
components and their interactions are termed
Earth System Models (ESM).

The high complexity of the Earth system not
only makes the approximation of the components
and their interaction very complicated and
involved (these code frameworks are huge with several persons' effort put in).
The actual computational power necessary to run these simulations
is extreme (Extreme scale Computing / Exascale
Computing).

@@colbox-blue
**Remark 1:** It is important to understand that, even
today, not all processes are fully understood. Hence,
their effect needs to be modeled as good
as possible.
@@

@@colbox-blue
**Remark 2:** Sometimes, the physics is well understood.
However in many cases nature is "multi-scale".
This means, for instance, that the spectra of
spatial sales range from centimeter and meter
up to 100+ km.
We will learn that "resolution"
makes simulations expensive to run on a computer.
The resolution necessary to resolve all scales
is unfortunately prohibitive, even on the
most powerful super computers. Hence, all
simulations are under-resolved. Hence,
effects from small un-resolved scales are
missing. It is, therefore, necessary to model
the missing so called subgrid scales.
@@

@@colbox-blue
**Remark 3:** In the climate weather prediction
community, the modeling of components, processes
and subgrid sale effects is often called "parametrization".
The "model" often refers to the whole package:
PDEs+Numerics+Implementation (code/software).
@@

As mentioned above, the good news is that depending
on the science question, some components processes
are more important than others. Hence, with
enough expertise it is feasible to choose a
subset of physics and consider simplified
climate models.


It thus is not surprising that there is a huge
collection of simplified climate models developed
by different researchers in the last decades.
Here is, for instance, a figure that lists a matrix
of simplifications with focus on Atmosphere+Hydrosphere:

![](/assets/milestone1/Models.png)
* Figure from lecture notes: [Stocker, "Introduction to Climate Modeling". Universität Bern](https://climatehomes.unibe.ch/~stocker/papers/stocker18icm.pdf), published by [Springer](https://link.springer.com/book/10.1007/978-3-642-00773-6).

## How did we choose a model for this course?

We used the following criteria to find a model for this course:

* It must be feasible to implement the model from scratch in one semester by students (ESM is out).
* We want to have a 2D grid to include earth surface modeling aspects that depend on the geography. Moreover, the topic of how to mesh the spheres is an important decision for a grid based model in GCMs/ESMs.
* The physics, the numerics, and the computational aspects should be accessible to a broad range of students (math, physics, meteorology, geophysics, etc.)

**What model did we choose?**

From the criteria mentioned above, we
decided to go for a 2D Energy Balance Model (EBM),
which is a heavily simplified climate model
that estimates the average temperature
of the atmosphere at the surface of Earth using a second order _partial differential equation_.

In our course, we will mainly follow the paper (although we changed the numerics a little bit):
> [Zhuang, K., North, G. R., & Stevens, M. J. (2017). A NetCDF version of the two-dimensional energy balance model based on the full multigrid algorithm. SoftwareX, 6, 198-202.](https://www.sciencedirect.com/science/article/pii/S2352711017300262)


