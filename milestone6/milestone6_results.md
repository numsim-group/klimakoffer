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

\fig{/assets/milestone6/annual_temperature.gif}

~~~
<iframe src="/assets/milestone6/temperature.html" seamless width=800 height=600></iframe>
~~~
### Mean Temperature Plot using the 2D-EBM

Plot of the mean and average temperature calculated by using the 2D-EBM and then averaging over the results:

\input{plot:1}{/assets/scripts/milestone6.jl}

### Cologne Temperature using the 2D-EBM

Cologne temperature calculated using the 2D-EBM model:

\input{plot:2}{/assets/scripts/milestone6.jl}

### Temperature for the NASA CO2 Data

Annual Temperatures from 1959 to 2020 based on NASA CO2 data:

\input{plot:3}{/assets/scripts/milestone6.jl}

### Mean Temperature Plot for Ziegler's data

Plot of the mean and average temperature calculated by using the 2D-EBM with Ziegler's parameters and then averaging over the results:

\input{plot:4}{/assets/scripts/milestone6.jl}


## Files Download

The Python implementation of milestone 6 can be downloaded here (Julia version will be added soon):

* [Julia solution](/assets/scripts/milestone6.jl)
* [Python solution](/assets/scripts/milestone6.py)

## Scripts for Milestone 6

You can also check out our Python implementation of milestone 6 in this site.

### Julia implementation of milestone 6

\input{julia}{/assets/scripts/milestone6.jl} 

### Python implementation of milestone 6

\input{python}{/assets/scripts/milestone6.py}
