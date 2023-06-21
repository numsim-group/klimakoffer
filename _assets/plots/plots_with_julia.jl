using DelimitedFiles
using Plots
using LaTeXStrings


# Data digitalised with https://apps.automeris.io/wpd/

alb_feb_avg_sky = readdlm("feb_avg.csv")
alb_feb_clear_sky = readdlm("feb_clear.csv")

figsize = (600,640) # Width x Height

# Creating plot 1 for average and clear sky albedo in February

p1 = plot(alb_feb_avg_sky[:,1], alb_feb_avg_sky[:,2], 
        xlabel = "Latitude [°]", ylabel = "Albedo [%]",
        xlims = [-70,70 ], ylims=[0,100],
        title = "Albedo vs Latitude (February)", label = "Average sky",
        size=figsize)

for i in [alb_feb_clear_sky] # Easy to extend the plotted data; just add it to the list and create a string list as well for the label
    plot!(p1,i[:,1],i[:,2], label = "Clear sky")
end

savefig(p1,"AlbedoLatitude.png")

# Creating plot comparing the modelling function with the data

x = LinRange(-70,70,140)
y = 100 .*(0.29 .+ (0.12*0.5) .* (3 .* sin.(deg2rad.(x)).^2 .- 1)) # Modelling function

p2 = plot(alb_feb_avg_sky[:,1], alb_feb_avg_sky[:,2], 
    xlabel = "Latitude [°]", ylabel = "Albedo [%]",
    xlims = [-70,70 ], ylims=[0,100],
    title = "Modelling albedo between the pole regions", label = "Average sky (February)",
    size = figsize)

plot!(p2, x,y,label = "0.3 + 0.12 p(θ)", color = "green")

savefig(p2,"ModellingAlbedo.png")



# Creating Radiative forcing plots

figsize_rad_feed = (740,670)


co2_ipcc = readdlm("co2_forcing_ipcc.csv")

p3 = plot(co2_ipcc[:,1], co2_ipcc[:,2], 
    xlabel = L"CO_2 \quad [ppm]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 1100], ylims=[0,10],
    title = L"Radiative\, forcing\, -\, CO_2", label = "",
    size = figsize_rad_feed
    )
savefig(p3, "CO2_forcing.png")



ch4_ipcc = readdlm("ch4_forcing_ipcc.csv")

p4 = plot(ch4_ipcc[:,1], ch4_ipcc[:,2], 
    xlabel = L"CH_4 \quad [ppb]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 5500], ylims=[0,2],
    title = L"Radiative\, forcing\, -\, CH_4", label = "",
    size = figsize_rad_feed
    )
savefig(p4, "CH4_forcing.png")


n2o_ipcc = readdlm("n2o_forcing_ipcc.csv")

p5 = plot(n2o_ipcc[:,1], n2o_ipcc[:,2], 
    xlabel = L"N_2\,O \quad [ppb]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [200, 600], ylims=[0,1],
    title = L"Radiative\, forcing\, -\, N_2\,O", label = "",
    size = figsize_rad_feed
    )
savefig(p5, "N2O_forcing.png")