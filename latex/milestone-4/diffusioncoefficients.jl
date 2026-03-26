using DelimitedFiles

"""
calc_diffusion_coefficients()
Calculate the diffusion coefficients with constants:
coeff_eq     = 0.65  coefficinet for diffusion at equator
coeff_ocean  = 0.40  coefficient for ocean diffusion
coeff_land   = 0.65  coefficinet for land diffusion
coeff_landNP = 0.28  coefficinet for land diffusion (north pole)
coeff_landSP = 0.20  coefficinet for land diffusion (south pole)
"""
function calc_diffusion_coefficients(geography,nlongitude=128,nlatitude=65)

  coeff_eq     = 0.65 
  coeff_ocean  = 0.40 
  coeff_land   = 0.65 
  coeff_landNP = 0.28 
  coeff_landSP = 0.20 

  diffusion = zeros(Float64,nlongitude,nlatitude)

  j_equator = div(nlatitude,2) + 1

  for j = 1:nlatitude
      theta = pi*real(j-1)/real(nlatitude-1)
      colat = sin(theta)^5

      for i = 1:nlongitude
          let geo = geography[i,j]
              if geo == 5 # oceans
                  diffusion[i,j] = (coeff_eq-coeff_ocean)*colat + coeff_ocean
              else # land, sea ice, etc
                  if j <= j_equator # northern hemisphere
                      diffusion[i,j] = (coeff_land-coeff_landNP)*colat + coeff_landNP
                  else # southern hemisphere
                      diffusion[i,j] = (coeff_land-coeff_landSP)*colat + coeff_landSP
                  end
                end
          end
      end
  end

  return diffusion
end


function milestone4()
    include(joinpath(@__DIR__,"..","milestone-1","geo.jl"))

    nlongitude = 128
    nlatitude = 65 

    geography = read_geography(joinpath(@__DIR__,"input","The_World128x65.dat"),nlongitude,nlatitude)

    return calc_diffusion_coefficients(geography,nlongitude,nlatitude)
end
