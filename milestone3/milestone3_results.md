+++
title = "Milestone 3"
hascode = true
rss = "Description"
rss_title = "Milestone 3"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++


# Milestone 3 - Project Results

\toc

## Expected Result

This is an example for the expected result:

### Solar Forcing vector
Spatial averages for every time step as a vector in a .txt file:\\
[(Download *solar\_forcing\_averages.txt*)](/assets/milestone3/solar_forcing_averages.txt)

### Annual Temperature Plot
Calculated with the forward Euler method. Plot for the backward Euler is almost identical.

\input{plot}{/assets/scripts/milestone3.jl}
