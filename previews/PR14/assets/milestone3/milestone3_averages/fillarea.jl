# This file was generated, do not modify it. # hide
for j=2:n_latitude-1
    theta_j = delta_theta * (j - 1)
    area[j] = sin(0.5 * delta_theta) * sin(theta_j) / n_longitude
end
# We print the area array to check that everything is right
println(area)