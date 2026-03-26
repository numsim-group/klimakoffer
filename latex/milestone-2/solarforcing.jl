using DelimitedFiles
using Interpolations
using Plots
pyplot()


"""
calc_albedo()
Calculates the albedo depending on the geography as
proposed by Zhuang et al. using a Legendre polynomial of second order
"""
function calc_albedo(geography,nlongitude=128,nlatitude=65)
    result = zeros(Float64,nlongitude,nlatitude)
    dtheta= pi/(nlatitude-1.0)
  
    for j in 1:nlatitude
      theta = 0.5*pi-dtheta*(j-1)
      sintheta = sin(theta)
      legendrepol = (3*sintheta^2-1)*0.5
      for i in 1:nlongitude
        let geo = geography[i,j]
          if geo == 1                           # land
            result[i,j] = 0.3+0.12*legendrepol
          elseif geo == 2                       # sea ice
            result[i,j] = 0.6
          elseif geo == 3                       # snow cover
            result[i,j] = 0.75
          elseif geo in (5,6,7,8)               # oceans
            result[i,j] = 0.29+0.12*legendrepol
          end
        end
      end
    end
    return result
  end
  

  """
  calc_heat_capacity()
  Calculates the heat capacity depending on the geography  
  """
function calc_heat_capacity(geography)

  # Depths 
  depth_mixed_layer = 70.0    # meters
  depth_soil = 2.0            # meters
  depth_seaice = 2.5          # meters
  depth_snow = 2.0            # meters
  layer_depth = 0.5           # kilometers

  # Physical properties of atmosphere
  rho_atmos = 1.293           # kg m^-3  dry air (STP)
  csp_atmos = 1005.0          # J kg^-1 K^-1 (STP)
  scale_height = 7.6          # kilometers
      
  # Physical properties of water
  rho_water = 1000.0          # kg m^-3
  csp_water = 4186.0          # J kg^-1 K^-1

  # Physical properties of soil 
  rho_soil = 1100.0           # kg m^-3   
  csp_soil = 850.0            # J kg^-1 K^-1
    
  # Physical properties of sea ice
  rho_sea_ice = 917.0         # kg m^-3
  csp_sea_ice = 2106.0        # J kg^-1 K^-1  

  # Physical properties of snow covered surface
  rho_snow = 400.0            # kg m^-3
  csp_snow = 1900.0           # J kg^-1 K^-1
      
  # Other constants  
  sec_per_yr = 3.15576e7      # seconds per year


  # atmosphere with exponentially decaying density
  sum = 0.0
  for n in 1:10
    z = (0.25 + layer_depth*real(n-1))/scale_height
    sum = sum + exp(-z)
  end

  c_atmos  	= csp_atmos*layer_depth*1000.0*rho_atmos*sum/sec_per_yr
  c_soil   	= depth_soil*rho_soil*csp_soil/sec_per_yr 
  c_seaice 	= depth_seaice*rho_sea_ice*csp_sea_ice/sec_per_yr
  c_snow   	= depth_snow * rho_snow * csp_snow/sec_per_yr
  c_mixed_layer = depth_mixed_layer*rho_water*csp_water/sec_per_yr

  # define heatcap
  heatcap = zeros(size(geography,1),size(geography,2))

  # Assign the correct value of the heat capacity of the columns
  for j in 1:size(geography,2)
    for i in 1:size(geography,1)
      geo  = geography[i,j]
      if geo == 1                            # land
        heatcap[i,j] = c_soil + c_atmos  
      elseif geo == 2                        # perennial sea ice
        heatcap[i,j] = c_seaice + c_atmos
      elseif geo == 3                        # permanent snow cover 
        heatcap[i,j] = c_snow + c_atmos         
      elseif geo == 5                        # ocean
        heatcap[i,j] = c_mixed_layer + c_atmos
      end                           
    end
  end  

  return heatcap
end


"""
calc_solar_forcing()
* Default s0 is 1371.685 [W/m²] (Current solar constant)
* Default orbital parameters of correspond to year 1950 AD:
    ecc = 0.016740             
    ob  = 0.409253
    per = 1.783037  
"""
function calc_solar_forcing(co_albedo, nlongitude, nlatitude, ntimesteps, yr=0; solar_cycle=false, s0=1371.685, orbital=false, ecc=0.016740, ob=0.409253, per=1.783037)
  # Calculate the sin, cos, and tan of the latitudes of Earth from the
  # colatitudes, calculate the insolation

  dy = pi/(nlatitude-1.0)

  siny = zeros(Float64,nlatitude)
  cosy = zeros(Float64,nlatitude)
  tany = zeros(Float64,nlatitude)

  for j in 1:nlatitude
    lat = pi/2.0 - dy*(j-1)    # latitude in radians
    siny[j] = sin(lat)
    if j == 1
      cosy[j] = 0.0
      tany[j] = 1000.0
    elseif j == nlatitude
      cosy[j] = 0.0
      tany[j] = -1000.0
    else
      cosy[j] = cos(lat)
      tany[j] = tan(lat)
    end
  end
  lambda = read_lambda(joinpath(@__DIR__,"input","True_Longitude.dat"))
  solar = _calc_insolation(ob, ecc, per, ntimesteps, nlatitude, siny, cosy, tany, s0, lambda)
  
  solar_forcing = zeros(Float64,nlongitude,nlatitude,ntimesteps)
  # calcualte the seasonal forcing
  for ts in 1:ntimesteps
    for j in 1:nlatitude
      for i in 1:nlongitude
        solar_forcing[i,j,ts] = solar[j,ts]*co_albedo[i,j]
      end
    end
  end
  return solar_forcing

end


"""
_calc_insolation()
auxiliar function of calc_solar_forcing
"""
function _calc_insolation(ob, ecc, per, nt, nlat, siny, cosy, tany, s0, lambda)

  eccfac = 1.0 - ecc^2
 
  solar = zeros(Float64,nlat,nt)


  #  Compute the average daily solar irradiance as a function of
  #  latitude and longitude (time)

  for n in 1:nt
    nu = lambda[n] - per
    rhofac = ((1.0- ecc*cos(nu))/eccfac)^2
    sindec = sin(ob)*sin(lambda[n])
    cosdec = sqrt(1.0-sindec^2)
    tandec = sindec/cosdec
    for j in 1:nlat
      z = -tany[j]*tandec
      if z >= 1.0    # polar latitudes when there is no sunrise (winter)
        solar[j,n] = 0.0
      else
        if z <= -1.0                # when there is no sunset (summer)
          solar[j,n] = rhofac * s0 * siny[j] * sindec
        else
          h_zero = acos(z)
          solar[j,n] = rhofac/pi*s0*(h_zero*siny[j]*sindec+cosy[j]*cosdec*sin(h_zero))
        end
      end
    end
  end

  return solar
end


"""
read_lambda()
Read input data to get the true longitude lambda.
"""
function read_lambda(filepath=joinpath(@__DIR__,"input","True_Longitude.dat"))

    result = readdlm(filepath, Float64)
return result
end

function milestone2()

  include(joinpath(@__DIR__,"..","milestone-1","geo.jl"))

  nlongitude = 128
  nlatitude = Int(nlongitude/2+1) # 65

  ntimesteps = 48

  geography = read_geography(joinpath(@__DIR__,"input","The_World128x65.dat"),nlongitude,nlatitude)
  albedo = calc_albedo(geography,nlongitude,nlatitude)
  co_albedo = 1.0 .- albedo

  heatcapacity = calc_heat_capacity(geography)
  solar_forcing = calc_solar_forcing(co_albedo, nlongitude, nlatitude, ntimesteps)


  # Read in the outline of the continents.
  outline = read_geography(joinpath(@__DIR__,"input","The_World_Outline128x65.dat"),128,65);
  xol,yol,zol = apply_robinson_projection(outline);

  x,y,z = apply_robinson_projection(solar_forcing[:,:,1]);

  anim = @animate for t in 1:size(solar_forcing,3)

      title = "Solar Forcing of the Earth: t = " * lpad(string(Int(floor(t*365/48))), 6, ' ') * " days";

      data = solar_forcing[:,:,t];
      z = reverse(reverse(data),dims=1);
      contourf(x,y,z,
          clims=(0,400),
          levels=LinRange(0,400,20),
          aspect_ratio=1.0,
          title=title,
          xlabel="longitude [°]",
          ylabel="latitude [°]",
          colorbar_title="scale of solar forcing",
          size=(1024,640));
      
      contour!(xol,yol,zol,c=:grays);
      
  end;

  gif(anim, "solar_forcing.gif", fps = 7)
end
