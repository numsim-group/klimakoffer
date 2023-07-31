+++
title = "Milestone 4"
hascode = true
rss = "Description"
rss_title = "Milestone 4"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 4 - The pointwise zero-dimensional Energy Balance Model

## The pointwise zero-dimensional Energy Balance Model

[(Download description in PDF format)](/assets/milestone4/description.pdf)\\
[(Download *The_World128x65.dat*)](/assets/scripts/input/The_World128x65.dat)\\

In the previous milestone you have implemented the zero-dimensional EBM for the global mean temperature. In this milestone you will use this model to compute the mean temperature for the northern and southern hemisphere. Furthermore you will implement the EBM model pointwise for each grid point.

1. Implement the functions *calc\_mean\_north* and *calc\_mean\_south* analogous to the function *calc\_mean* from milestone 3 to calculate the mean value on the southern and northern hemisphere of a given parameter for an input matrix of values of that paramer in each grid point.

2. Use the *calc\_equilibrium* function from milestone 3 to compute the EBM for the northern and southern hemisphere by substituting the *calc\_mean* function with the two functions from task 1. What differences do you notice to the mean temperature result from milestone 3?

3. Now use the *calc\_equilibrium* function to compute the EBM for each grid point by looping over every grid point and substituting the mean values with the values on each grid point. Plot the results as an animation similar to milestone 2.

4. Use the results of the pointwise calculation in task 3 to calculate the mean temperatures over the whole globe, the northern hemisphere and the southern hemisphere using the functions *calc\_mean*, *calc\_mean\_north* and *calc\_mean\_south*. What differences are there to the results of the EBM calculations that are based on the mean values from the start?

5. What is the temperature in Cologne according to the pointwise EBM? Compare it to the other results. What is noticeable?
