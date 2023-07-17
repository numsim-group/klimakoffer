using DelimitedFiles
using Plots
using LaTeXStrings

linewidth = 2 # linewidth of plots
markersize = 10 # Radius of scatter plot symbols in px


# Data digitalised with https://apps.automeris.io/wpd/

alb_feb_avg_sky = readdlm("feb_avg.csv")
alb_feb_clear_sky = readdlm("feb_clear.csv")

figsize = (600,640) # Width x Height

# Correct x axis considering that the original plot has a sinusoidal axis... And the digitalization was done using reference x values of -60 and 60!
m = sin(deg2rad(60)) / 60
for i in 1:size(alb_feb_avg_sky)[1]
    alb_feb_avg_sky[i,1] *= m
end
for i in 1:size(alb_feb_clear_sky)[1]
    alb_feb_clear_sky[i,1] *= m
end

# Creating plot 1 for average and clear sky albedo in February
p1 = plot(alb_feb_avg_sky[:,1], alb_feb_avg_sky[:,2], 
        xlabel = "Latitude [°]", ylabel = "Albedo [%]",
        ylims=[0,100], #xlims = [-70,70 ], 
        title = "Albedo vs Latitude (February)", label = "Average sky",
        size=figsize,
        linewidth = linewidth)

plot!(p1,alb_feb_clear_sky[:,1],alb_feb_clear_sky[:,2], label = "Clear sky", linewidth = linewidth)

# Draw ticks that are consistent with the axis in sinusoidal scale, but write the angle in degrees!
plot!(xticks = (sin.(deg2rad.([-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])), ["-90", "", "-60", "", "-30", "", "0", "", "30", "", "60", "", "90"]))

savefig(p1,"AlbedoLatitude.png")

# Creating plot comparing the modelling function with the data

x = LinRange(-1,1,140)
y = 100 .*(0.29 .+ (0.12*0.5) .* (3 .* x.^2 .- 1)) # Modelling function


p2 = plot(alb_feb_avg_sky[:,1], alb_feb_avg_sky[:,2], 
    xlabel = "Latitude [°]", ylabel = "Albedo [%]",
    ylims=[0,100],
    title = "Modelling albedo between the pole regions", label = "Average sky (February)",
    size = figsize,
    linewidth = linewidth)

plot!(p2, x,y,label = "0.29 + 0.12 p(θ)", color = "green", linewidth = linewidth)
plot!(xticks = (sin.(deg2rad.([-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])), ["-90", "", "-60", "", "-30", "", "0", "", "30", "", "60", "", "90"]))

savefig(p2,"ModellingAlbedo.png")



# Creating Radiative forcing plots

figsize_rad_feed = (740,670)
fit_label = L"Fit\; to\; NBM\; data"
data_label = L"NBM\; (Narrow \;Band\; Model)\; results"

co2_ipcc = readdlm("co2_forcing_ipcc.csv")
co2_data = readdlm("co2_forcing_data.csv")


p3 = plot(co2_data[:,1],co2_data[:,2], seriestype=:scatter, shape = :+, 
    xlabel = L"CO_2 \quad [ppmv]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 1100], ylims=[0,10],
    title = L"Radiative\, forcing\, -\, CO_2", 
    label = data_label,
    size = figsize_rad_feed,
    markersize = markersize
    )

plot!(p3,co2_ipcc[:,1], co2_ipcc[:,2],label = fit_label, linewidth = linewidth)

savefig(p3, "CO2_forcing.png")



ch4_ipcc = readdlm("ch4_forcing_ipcc.csv")
ch4_data = readdlm("ch4_forcing_data.csv")

p4 = plot(ch4_data[:,1], ch4_data[:,2],seriestype=:scatter, shape = :+,
    xlabel = L"CH_4 \quad [ppbv]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 5500], ylims=[0,2],
    title = L"Radiative\, forcing\, -\, CH_4", label=data_label,
    size = figsize_rad_feed,
    markersize = markersize
    )

plot!(p4,ch4_ipcc[:,1], ch4_ipcc[:,2],label = fit_label, linewidth = linewidth)
savefig(p4, "CH4_forcing.png")


n2o_ipcc = readdlm("n2o_forcing_ipcc.csv")
n2o_data = readdlm("n2o_forcing_data.csv")

p5 = plot(n2o_data[:,1], n2o_data[:,2], seriestype=:scatter, shape=:+,
    xlabel = L"N_2\,O \quad [ppbv]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 600], ylims=[0,1],
    title = L"Radiative\, forcing\, -\, N_2\,O", label = data_label,
    size = figsize_rad_feed,
    markersize = markersize
    )

plot!(p5,n2o_ipcc[:,1],n2o_ipcc[:,2], label = fit_label, linewidth = linewidth)
savefig(p5, "N2O_forcing.png")


# Budyko's Linear Model

textbook_curve = readdlm("radiation_physical_model.csv")
fitting_curve = readdlm("radiation_fitting_curve.csv")
points = readdlm("radiation_points.csv")

p6 = plot(
    textbook_curve[:,1],
    textbook_curve[:,2],
    label = L"\sigma T^4",
    xlabel = L"Surface \;  Temperature \;  [°C]",
    ylabel = L"OLR \;[W/m^2]",
    title = L"Outgoing\; Longwave\; Radiation\; vs.\; Temperature",
    size = (800,600),
    xticks = LinRange(-40,40,9),
    legendfontsize = 12,
    titlefontsize = 20,
    guidefontsize = 15,
    linewidth = linewidth
)

plot!(
    p6,
    points[:,1],
    points[:,2],
    seriestype = :scatter,
    shape = :+,
    label = "",
    markersize = 5 # We manually set the markersize here, so the plot does not look super crammed. 
)

plot!(p6,
    fitting_curve[:,1],
    fitting_curve[:,2],
    label = L"202.1 + 1.90 T",
    linewidth = linewidth
)

savefig(p6, "OutgoingLongwaveRadiation.png")



# Heat Transfer Plot

heat_ocean = readdlm("heat_transfer_ocean.csv")
heat_total = readdlm("heat_transfer_total.csv")
heat_atmos = readdlm("heat_transfer_atmosphere.csv")
p7 = plot(
    heat_total[:,1],
    heat_total[:,2],
    title = "Heat transport in the Northern Hemisphere",
    xlabel = "Latitude [°]",
    ylabel = "Northward heat transport [10³ TW]",
    label = "Total",
    size = (600,400),
    xticks = LinRange(0,90,10),
    yticks = [-1, 0, 1, 2, 3, 4, 5, 6],
    linewidth = linewidth
)

plot!(
    p7,
    heat_ocean[:,1],
    heat_ocean[:,2],
    label = "Ocean",
    linewidth = linewidth
)

plot!(
    p7,
    heat_atmos[:,1],
    heat_atmos[:,2],
    label = "Atmosphere",
    linewidth = linewidth
)

savefig(p7, "heat_transfer_north.png")