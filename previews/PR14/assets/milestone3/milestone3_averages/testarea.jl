# This file was generated, do not modify it. # hide
field = ones(Float64,n_latitude,n_longitude)
println("Area : ", calc_mean(field, area, n_latitude, n_longitude))
println("Error: ", calc_mean(field, area, n_latitude, n_longitude)-1)