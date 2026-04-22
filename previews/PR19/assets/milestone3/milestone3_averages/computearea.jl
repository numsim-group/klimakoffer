# This file was generated, do not modify it. # hide
function calc_mean(field, area, n_latitude, n_longitude)
  if size(field)[1] !== n_latitude || size(field)[2] !== n_longitude
    error("field and area sizes do not match")
  end

  # Initialize mean with the values at the poles
  mean = area[1]*field[1,1] + area[end]*field[end,end]

  for j in 2:n_latitude-1
      for i in 1:n_longitude
          mean += area[j] * field[j,i]
      end 
  end

  return mean
end