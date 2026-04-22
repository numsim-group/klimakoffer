# This file was generated, do not modify it. # hide
delta_theta = pi / (n_latitude - 1)
area[1] = 0.5 * (1 - cos(0.5 * delta_theta))
area[n_latitude] = area[1]