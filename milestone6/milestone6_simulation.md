+++
title = "Milestone 6"
hascode = true
rss = "Description"
rss_title = "Milestone 6"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 6 - Climate Simulations

\toc

## Standard use cases with Klimakoffer

The general idea of the 2D EBM and its numerical realization, Klimakoffer, is to simulate the annual equilibrium temperature for a given value of the CO2 concentration. For instance, one could directly simulate the surface temperatures in the year 1950. 

Due to the relatively general and complex solar forcing model, it is also possible to adapt the orbital parameters to pre-historic times, e.g., 9ka BP or 21ka BP (kilo years before present) as demonstrated in our main reference 
> [Zhuang, K., North, G. R., & Stevens, M. J. (2017). A NetCDF version of the two-dimensional energy balance model based on the full multigrid algorithm. SoftwareX, 6, 198-202.](https://www.sciencedirect.com/science/article/pii/S2352711017300262)

The main use case is to investigate the CO2 response of the model, i.e., what happens when the CO2 value is doubled. How does the temperature field change, how does the average Earth temperature change. 


## Further capabilities of Klimakoffer

An immediate application with Klimakoffer is to investigate the sensitivity of the surface temperature field for changing geography. We stress, that most crucial parameters directly depend on the geography. Hence, changing the geography will cascade through the model and adapt all crucial parameters accordingly. Thus, another capability of Klimakoffer is to investigate the impact of changing geography. We stress, however, that the current model is not capable of a (non-linear) feedback mechanism, i.e., changing geography throughout the year depending on the temperature. This would require some (slight) modifications of the numerics as the Jacobian matrix would change in time and depend explicitly on the temperature field. In contrast, changing the fixed geography in the beginning of the simulation is an option and would already give some insight on the impact of, for instance, ice and snow distribution on Earth. 

Another straightforward application is to investigate the temperature response of Klimakoffer throughout the years with changing (increasing) CO2 concentration. Data of average CO2 increase is available, e.g., from NASA. 

\fig{/assets/milestone6/co2_nasa.png}
* **Source:** NASA ([https://climate.nasa.gov/vital-signs/carbon-dioxide/](https://climate.nasa.gov/vital-signs/carbon-dioxide/))

One would use the Klimakoffer in a sequential way, i.e., a chain of equilibrium simulations where for each year and given value of CO2 the equilibrium in the year is found and chained together. 

Last but not least, we put a lot of emphasis during the lecture and the modeling sessions to stress that all parameters are empirical and that several authors choose different fitting methods, resulting in another choice of parameters. Thus, another capability of the current Klimakoffer is to investigate the sensitivity of the resulting temperature fields to the choice of certain parameters in the model.