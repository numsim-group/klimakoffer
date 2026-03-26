using DelimitedFiles

"""
calc_area()
Calculates fractional area for the cells depending on latitude
"""
function calc_area(nx=128)

    ny = Int(nx/2+1)

    h = pi / (ny-1)          # Uniform grid size in radians


    area = zeros(Float64,ny) 

    for j=2:ny-1
        area[j] =  sin(0.5*h)*sin(h*(j-1))/nx
    end

    # Fractional area for the poles
    area[1]  = 0.5*(1 - cos(0.5*h))
    area[ny] = area[1]

    return area
end

"""
calc_mean_quantities()
Calculates average values for alebdo, heat capacity and solar forcing
"""
function calc_mean_quantities(nlongitude=128,nlatitude=65,ntimesteps=48)
    geography = read_geography(joinpath(@__DIR__,"input","The_World128x65.dat"),nlongitude,nlatitude)
    albedo = calc_albedo(geography,nlongitude,nlatitude)
    heat_capacity = calc_heat_capacity(geography)
    area = calc_area(nlongitude)
    co_albedo = 1.0 .- albedo
    solar_forcing = calc_solar_forcing(co_albedo, nlongitude, nlatitude, ntimesteps)

    # Compute mean solar forcing at each time step
    mean_solar_forcing = zeros(Float64,ntimesteps)
    for t in 1:ntimesteps
        mean_solar_forcing[t] = calc_mean(solar_forcing[:,:,t],area)
    end 

    return calc_mean(albedo,area), calc_mean(heat_capacity,area), mean_solar_forcing
end 

"""
calc_mean(matrix, area)
  Calculates average values for a matrix of values on the globe grid and an area array
  * This routine uses the first entry of the matrix for the north pole and the last entry for the south pole
"""
function calc_mean(matrix, area)
  # Get size of matrix
  nlongitude = size(matrix,1)
  nlatitude = size(matrix,2)

  if size(area)[1] !== nlatitude
    error("Matrix and area sizes do not match")
  end

  # Initialize mean with the values at the poles
  mean = area[1]*matrix[1,1] + area[end]*matrix[end,end]

  for j in 2:nlatitude-1
      for i in 1:nlongitude
          mean += area[j] * matrix[i,j]
      end 
  end

  return mean
end 

"""
    compute_equilibrium!(...)
* rel_error is the tolerance for global temperature equilibrium (default is 2e-5).
* max_years is the maximum number of annual cycles to be computed when searching for equilibrium
"""
function compute_equilibrium!(nlongitude=128, nlatitude=65, num_steps_year=48, max_years=100, rel_error=2e-5, verbose=true)

    annual_temperature = fill(5.0,num_steps_year)

    albedo, heat_capacity, solar_forcing = calc_mean_quantities(nlongitude, nlatitude, num_steps_year)
   
    radiative_cooling_co2= calc_radiative_cooling_co2(315.0)          #A
    radiative_cooling_feedback = 2.15                                 #B

    average_temperature = average_temperature_old = annual_temperature[num_steps_year]

    if verbose
      println("year","  ","Average Temperature")
      println(0,"  ",average_temperature_old)
    end

    for year in 1:max_years
        average_temperature = 0.0
        for time_step in 1:num_steps_year
            old_time_step = (time_step == 1) ? num_steps_year : time_step - 1
            annual_temperature[time_step]= (0.5*(solar_forcing[time_step]+solar_forcing[old_time_step])/num_steps_year-radiative_cooling_co2/num_steps_year+(heat_capacity-0.5*radiative_cooling_feedback/num_steps_year)*annual_temperature[old_time_step])/(heat_capacity+0.5*radiative_cooling_feedback/num_steps_year)
          
            average_temperature += annual_temperature[time_step]
        end

        average_temperature = average_temperature/num_steps_year

        if verbose
          println(year,"  ",average_temperature)
        end
        
        if (abs(average_temperature-average_temperature_old)<rel_error)
            if verbose
              println("EQUILIBRIUM REACHED!")
            end
            break
        end
        
        average_temperature_old = average_temperature
    end

    return average_temperature
end



"""
calc_radiative_cooling_co2()
Computes the CO2 parameter depending on the CO2 concentration
* Default CO2 concentration is 315 ppm (equivalent to year 1950)
"""
function calc_radiative_cooling_co2(co2_concentration=315.0)

  co2_concentration_base = 315.0
  radiative_cooling_co2_base = 210.3

  radiative_cooling_co2=radiative_cooling_co2_base-5.35*log(co2_concentration/co2_concentration_base)

  return radiative_cooling_co2
end


function milestone3()

  include(joinpath(@__DIR__,"..","milestone-1","geo.jl"))
  include(joinpath(@__DIR__,"..","milestone-2","solarforcing.jl"))

  nlongitude = 128
  nlatitude = 65 
  ntimesteps = 48

  return compute_equilibrium!(nlongitude,nlatitude,ntimesteps)
end
