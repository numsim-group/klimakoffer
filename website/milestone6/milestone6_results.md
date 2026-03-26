+++
title = "Milestone 6"
hascode = true
rss = "Description"
rss_title = "Milestone 6"
rss_pubdate = Date(2022, 5, 1)

tags = ["discrete ebm"]
+++


# Milestone 6 - Project Results

\toc

## Expected Results

This is an example for the expected results:

### Annual Temperature Animation

Annual Temperature calculated with the 2D-EBM model:

[Annual temperature animation](/assets/milestone_results/milestone6_annual_temperature.gif)

~~~
<iframe src="/assets/milestone6/temperature.html" seamless width=800 height=600></iframe>
~~~
### Mean Temperature Plot using the 2D-EBM

Plot of the mean and average temperature calculated by using the 2D-EBM and then averaging over the results:

[Mean temperature plot](/assets/milestone_results/milestone6_temperature.png)

### Cologne Temperature using the 2D-EBM

Cologne temperature calculated using the 2D-EBM model:

[Cologne temperature plot](/assets/milestone_results/milestone6_cologne.png)

### Temperature for the NASA CO2 Data

Annual Temperatures from 1959 to 2020 based on NASA CO2 data:

[Temperature for NASA CO2 data](/assets/milestone_results/milestone6_temperature_co2.png)

### Mean Temperature Plot for Ziegler's data

Plot of the mean and average temperature calculated by using the 2D-EBM with Ziegler's parameters and then averaging over the results:

\input{plot:4}{/assets/scripts/milestone6.jl}


## Files Download

The Python implementation of milestone 6 can be downloaded here (Julia version will be added soon):

* [Julia solution](/assets/scripts/milestone6.jl)
* [Python solution](/assets/scripts/milestone6.py)

## Scripts for Milestone 6

You can also check out our Julia and Python implementations of milestone 6 in this site.

### Julia implementation of milestone 6

\input{julia}{/assets/scripts/milestone6.jl}

### Python implementation of milestone 6

\input{python}{/assets/scripts/milestone6.py}
