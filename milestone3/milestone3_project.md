+++
title = "Milestone 3"
hascode = true
rss = "Description"
rss_title = "Milestone 3"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 3 - Solving a zero-dimensional Energy Balance Model

## Solving a zero-dimensional Energy Balance Model

[(Download description in PDF format)](/assets/milestone3/description.pdf)\\
[(Download *The_World128x65.dat*)](/assets/scripts/input/The_World128x65.dat)\\
[(Download *True_Longitude.dat*)](/assets/scripts/input/True_Longitude.dat)

In the previous milestones you have already implemented some routines that are part of an energy balance model. Now you have the opportunity to combine your results into a zero-dimensional ($0D$) energy balance model (EBM) and solve the model using an ODE solver. Note that for zero-dimensional EBMs, spatial dependence is neglected, and $T$ is a global time-dependent averaged surface temperature. Let the following equation be given for the energy balance model:

\begin{align}
\overline{C}\frac{\partial {T}_A(t)}{\partial t}+A+B{T}_A(t)=\overline{S_\text{sol}}(t),
\label{EBM}
\end{align}
where $\overline{C}$ is the mean heat capacity, ${T}_A$ is an approximation to the area-average of the time dependent temperature, $A$ the $CO_2$ dependent radiative cooling, $B = 2.15$ the feedback effects such as water vapor cycles, lapse rate and cloud cover, $\overline{S_\text{sol}}(t)$ is the area-averaged solar forcing.
To solve the zero-dimensional EBM you should proceed as described in the next steps.

1. Implement a function *calc_area*, which determines the area fraction of each grid cell as a function of latitude, with an equidistant step size $h=\frac{\pi}{64}$ for a $128 \times 65$ grid, and returns the values as a vector.
2. Implement a *calc_mean* function, using the results from the first task, which takes a matrix of values of a given parameter like albedo or solar forcing for each grid cell and returns the mean value over the whole area.
3. Write a function *calc\_radiative\_cooling\_co2*, which returns the $CO_2$ dependent radiative cooling $A$ for a given $CO_2$ concentration $c$. As default value use the $CO_2$ concentration of the year 1950: $c_0 = 315.0$. The radiative cooling $A$ is given by
\begin{align}
A(c) = 210.3-5.35 \log\left(\frac{c}{c_0}\right)
\label{A}
\end{align}
4. In a function *compute_equilibrium*, solve the above equation \eqref{EBM} using the forward Euler method and calculate an average temperature per year until equilibrium is reached.
5. Write a new function or adapt your existing one to add an option to solve the equation with the backward Euler method instead. Plot your results.
