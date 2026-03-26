using Plots
using LaTeXStrings
pythonplot()

lat = LinRange(-90,90,180)
lat_north = LinRange(0,90,90)
lat_south = LinRange(-90,0,90)
#print(lat)
ocean = 0.4 .+ (0.65-0.4) .* sin.(deg2rad.(lat.+90)).^5
north_nonocean = 0.28 .+ (0.65-0.28) .* sin.(deg2rad.(lat_north.+90)).^5
south_nonocean = 0.2 .+ (0.65-0.2) .* sin.(deg2rad.(lat_south.+90)).^5
nonocean = vcat(south_nonocean, north_nonocean)
plot(lat,ocean, levels=LinRange(-90, 90, 7), label="over ocean", xlabel = "latitude [°]", ylabel = L"\widetilde{D} \: \: [W/m^2/K]", xlims=[-90,90], ylims=[0,1], yticks = 6)
plot!(lat,nonocean, label="over land, snow, sea ice", linecolor="red")
