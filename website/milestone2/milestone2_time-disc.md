+++
title = "Milestone 2"
hascode = true
rss = "Description"
rss_title = "Milestone 2"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 2 - A note on the time discretization
\toc

In our simple climate model, we will use 48 time steps each year and set the first time step as the vernal equinox. For such a time-discretization with 48 time steps, the main astronomical events and their correspondence with time steps are listed in the following table:

|  Astronomical event |   Time step|
| ------------------ | -----------| 
|   Vernal Equinox   |   1        |  
|   Summer Solstice  |   13       | 
|   Autumnal Equinox |   25       | 
|   Winter Solstice  |   37       | 

In addition, the following table presents the time steps that correspond to each month:

| Month |     Time Steps  | 
|-------|-----------------|
| Jan   |  38, 39, 40, 41 |
| Feb   |  42, 43, 44, 45 |
| Mar   |  46, 47, 48   1 |
| Apr   |   2,  3,  4,  5 |
| May   |   6,  7,  8,  9 |
| Jun   |  10, 11, 12, 13 |
| Jul   |  14, 15, 16, 17 |
| Aug   |  18, 19, 20, 21 |
| Sep   |  22, 23, 24, 25 |
| Oct   |  26, 27, 28, 29 |
| Nov   |  30, 31, 32, 33 |
| Dec   |  34, 35, 36, 37 |
