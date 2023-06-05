+++
title = "Milestone 6"
hascode = true
rss = "Description"
rss_title = "Milestone 6"
rss_pubdate = Date(2022, 5, 1)

tags = ["linear systems", "linear solvers"]
+++


# Milestone 6 - Fixing Stability in the spatial 2D Energy Balance Model

## Fixing Stability in the spatial 2D Energy Balance Modell

[(Download description in PDF format)](/assets/milestone6/description.pdf)\\
[(Download *The\_World128x65.dat*)](/assets/milestone6/input/The_World128x65.dat)\\
[(Download *True\_Longitude.dat*)](/assets/milestone6/input/True_Longitude.dat)\\
[(Download *co2\_nasa.dat*)](/assets/milestone6/input/co2_nasa.dat)\\

1. In the previous milestone you have implemented the spatial operator of the diffusive EBM for two spatial dimensions. The investigations of the eigenvalues revealed that the stable time step of the explicit Euler method is prohibitively small. Hence, we focus on the implicit Euler method. Adapt the function *calc\_equilibrium\_2d* to use the backward Euler instead of the forward Euler time integration method. Visualize the results as an animation.

2. Use the functions *calc\_area\_mean*, *calc\_area\_mean\_north* and *calc\_area\_mean\_south* from previous milestones to compute the mean temperature values from the EBM in two space dimensions. Compare them with the results in milestone 4. Compare also the results for the temperature in Cologne. What is noticeable?

3. Write a function to compute equilibrium simulations based on NASA $CO_2$ data from 1980 to the present. Plot the results.

4. Now we want to investigate the sensitivity of the model to the parameters. Repeat the calculations of the temperature for the year 1950 with the parameters from Ziegler and Rehfeld 2021[^1]

\fig{/assets/milestone6/new_parameter.png}

Compare your results with the temperatures calculated in 1.

5. Compute the equilibrium temperature from 1980 to present with the NASA $CO_2$ data as in task 3 but with the new parameter. Compare your results.

[^1]: E. Ziegler and K.Rehfeld, TransEBM v. 1.0: description, tuning, and validation of a transient model of the Earth's energy balance in two dimensions, Geoscientific Model Development 14.5, 2021.