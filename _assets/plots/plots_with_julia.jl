using DelimitedFiles
using Plots
using LaTeXStrings


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
        size=figsize)

plot!(p1,alb_feb_clear_sky[:,1],alb_feb_clear_sky[:,2], label = "Clear sky")

# Draw ticks that are consistent with the axis in sinusoidal scale, but write the angle in degrees!
plot!(xticks = (sin.(deg2rad.([-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])), ["-90", "", "-60", "", "-30", "", "0", "", "30", "", "60", "", "90"]))

savefig(p1,"AlbedoLatitude.png")

# Creating plot comparing the modelling function with the data

x = LinRange(-70,70,140)
y = 100 .*(0.29 .+ (0.12*0.5) .* (3 .* deg2rad.(x).^2 .- 1)) # Modelling function

x =  m .* x

p2 = plot(alb_feb_avg_sky[:,1], alb_feb_avg_sky[:,2], 
    xlabel = "Latitude [°]", ylabel = "Albedo [%]",
    ylims=[0,100],
    title = "Modelling albedo between the pole regions", label = "Average sky (February)",
    size = figsize)

plot!(p2, x,y,label = "0.3 + 0.12 p(θ)", color = "green")
plot!(xticks = (sin.(deg2rad.([-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])), ["-90", "", "-60", "", "-30", "", "0", "", "30", "", "60", "", "90"]))

savefig(p2,"ModellingAlbedo.png")



# Creating Radiative forcing plots

figsize_rad_feed = (740,670)
fit_label = L"Fit\; to\; NBM\; data"
data_label = L"NBM\; (Narrow \;Band\; Model)\; results"

co2_ipcc = readdlm("co2_forcing_ipcc.csv")
co2_data = readdlm("co2_forcing_data.csv")


p3 = plot(co2_data[:,1],co2_data[:,2], seriestype=:scatter, shape = :+, 
    xlabel = L"CO_2 \quad [ppm]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 1100], ylims=[0,10],
    title = L"Radiative\, forcing\, -\, CO_2", 
    label = data_label,
    size = figsize_rad_feed
    )

plot!(p3,co2_ipcc[:,1], co2_ipcc[:,2],label = fit_label)

savefig(p3, "CO2_forcing.png")



ch4_ipcc = readdlm("ch4_forcing_ipcc.csv")
ch4_data = readdlm("ch4_forcing_data.csv")

p4 = plot(ch4_data[:,1], ch4_data[:,2],seriestype=:scatter, shape = :+,
    xlabel = L"CH_4 \quad [ppb]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 5500], ylims=[0,2],
    title = L"Radiative\, forcing\, -\, CH_4", label=data_label,
    size = figsize_rad_feed
    )

plot!(p4,ch4_ipcc[:,1], ch4_ipcc[:,2],label = fit_label)
savefig(p4, "CH4_forcing.png")


n2o_ipcc = readdlm("n2o_forcing_ipcc.csv")
n2o_data = readdlm("n2o_forcing_data.csv")

p5 = plot(n2o_data[:,1], n2o_data[:,2], seriestype=:scatter, shape=:+,
    xlabel = L"N_2\,O \quad [ppb]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 600], ylims=[0,1],
    title = L"Radiative\, forcing\, -\, N_2\,O", label = data_label,
    size = figsize_rad_feed
    )

plot!(p5,n2o_ipcc[:,1],n2o_ipcc[:,2], label = fit_label)
savefig(p5, "N2O_forcing.png")