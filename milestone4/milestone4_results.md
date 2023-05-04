+++
title = "Milestone 4"
hascode = true
rss = "Description"
rss_title = "Milestone 4"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++


# Milestone 4 - Project Results

\toc

## Expected Results

This is an example for the expected results:

### Mean Temperature Plot

Plot of the mean and average temperature calculated by averaging over the regions first and then using the 0D-EBM.

\input{plot:1}{/assets/scripts/milestone4.jl}

### Mean Temperature Plot using Pointwise EBM

Plot of the mean and average temperature calculated by using a 0D-EBM for every point of the discretization and then averaging over the results.

\input{plot:2}{/assets/scripts/milestone4.jl}

### Cologne Temperature

Cologne temperature calculated via the pointwise 0D-EBM

\input{plot:3}{/assets/scripts/milestone4.jl}

### Annual Temperature Animation
Annual temperature with the pointwise 0D-EBM:

\fig{/assets/milestone4/annual_temperature.gif}

