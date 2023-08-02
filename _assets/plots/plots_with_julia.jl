using DelimitedFiles
using Plots
using LaTeXStrings

linewidth = 2 # linewidth of plots
markersize = 10 # Radius of scatter plot symbols in px
guidefontsize = 14 # Size of text
dpi = 300 # resolution

# Data digitalised with https://apps.automeris.io/wpd/

alb_feb_avg_sky = readdlm(joinpath(@__DIR__,"feb_avg.csv"))
alb_feb_clear_sky = readdlm(joinpath(@__DIR__,"feb_clear.csv"))

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
        ylims=[0,100],
        # xlims = [-70,70 ], 
        title = "Albedo vs Latitude (February)", label = "Average sky",
        size=figsize,
        linewidth = linewidth,
        guidefontsize = guidefontsize,
        dpi = dpi)

plot!(p1,alb_feb_clear_sky[:,1],alb_feb_clear_sky[:,2], label = "Clear sky", linewidth = linewidth)

# Draw ticks that are consistent with the axis in sinusoidal scale, but write the angle in degrees!
plot!(xticks = (sin.(deg2rad.([-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])), ["-90", "", "-60", "", "-30", "", "0", "", "30", "", "60", "", "90"]))

savefig(p1,joinpath(@__DIR__,"AlbedoLatitude.png"))

# Creating plot comparing the modelling function with the data

x = LinRange(-1,1,140)
y = 100 .*(0.29 .+ (0.12*0.5) .* (3 .* x.^2 .- 1)) # Modelling function


p2 = plot(alb_feb_avg_sky[:,1], alb_feb_avg_sky[:,2], 
    xlabel = "Latitude [°]", ylabel = "Albedo [%]",
    ylims=[0,100],
    title = "Modelling albedo between the pole regions", label = "Average sky (February)",
    size = figsize,
    linewidth = linewidth,guidefontsize = guidefontsize, dpi = dpi)

plot!(p2, x,y,label = "0.29 + 0.12 p(θ)", color = "green", linewidth = linewidth)
plot!(xticks = (sin.(deg2rad.([-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])), ["-90", "", "-60", "", "-30", "", "0", "", "30", "", "60", "", "90"]))

savefig(p2,joinpath(@__DIR__,"ModellingAlbedo.png"))



# Creating Radiative forcing plots

figsize_rad_feed = (740,670)
fit_label = L"Fit\; to\; NBM\; data"
data_label = L"NBM\; (Narrow \;Band\; Model)\; results"

co2_ipcc = readdlm(joinpath(@__DIR__,"co2_forcing_ipcc.csv"))
co2_data = readdlm(joinpath(@__DIR__,"co2_forcing_data.csv"))


p3 = plot(co2_data[:,1],co2_data[:,2], seriestype=:scatter, shape = :+, 
    xlabel = L"CO_2 \quad [ppmv]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 1100], ylims=[0,10],
    title = L"Radiative\, forcing\, -\, CO_2", 
    label = data_label,
    size = figsize_rad_feed,
    markersize = markersize,guidefontsize = guidefontsize, dpi = dpi
    )

plot!(p3,co2_ipcc[:,1], co2_ipcc[:,2],label = fit_label, linewidth = linewidth)

savefig(p3, joinpath(@__DIR__,"CO2_forcing.png"))



ch4_ipcc = readdlm(joinpath(@__DIR__,"ch4_forcing_ipcc.csv"))
ch4_data = readdlm(joinpath(@__DIR__,"ch4_forcing_data.csv"))

p4 = plot(ch4_data[:,1], ch4_data[:,2],seriestype=:scatter, shape = :+,
    xlabel = L"CH_4 \quad [ppbv]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 5500], ylims=[0,2],
    title = L"Radiative\, forcing\, -\, CH_4", label=data_label,
    size = figsize_rad_feed,
    markersize = markersize,guidefontsize = guidefontsize, dpi = dpi
    )

plot!(p4,ch4_ipcc[:,1], ch4_ipcc[:,2],label = fit_label, linewidth = linewidth)
savefig(p4, joinpath(@__DIR__,"CH4_forcing.png"))


n2o_ipcc = readdlm(joinpath(@__DIR__,"n2o_forcing_ipcc.csv"))
n2o_data = readdlm(joinpath(@__DIR__,"n2o_forcing_data.csv"))

p5 = plot(n2o_data[:,1], n2o_data[:,2], seriestype=:scatter, shape=:+,
    xlabel = L"N_2\,O \quad [ppbv]", ylabel = L"Radiative\, forcing \quad [W/m^2]",
    xlims = [0, 600], ylims=[0,1],
    title = L"Radiative\, forcing\, -\, N_2\,O", label = data_label,
    size = figsize_rad_feed,
    markersize = markersize,guidefontsize = guidefontsize, dpi = dpi
    )

plot!(p5,n2o_ipcc[:,1],n2o_ipcc[:,2], label = fit_label, linewidth = linewidth)
savefig(p5, joinpath(@__DIR__,"N2O_forcing.png"))


# Budyko's Linear Model

textbook_curve = readdlm(joinpath(@__DIR__,"radiation_physical_model.csv"))
fitting_curve = readdlm(joinpath(@__DIR__,"radiation_fitting_curve.csv"))
points = readdlm(joinpath(@__DIR__,"radiation_points.csv"))

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
    guidefontsize = guidefontsize, dpi = dpi,
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

savefig(p6, joinpath(@__DIR__,"OutgoingLongwaveRadiation.png"))



# Heat Transfer Plot

heat_ocean = readdlm(joinpath(@__DIR__,"heat_transfer_ocean.csv"))
heat_total = readdlm(joinpath(@__DIR__,"heat_transfer_total.csv"))
heat_atmos = readdlm(joinpath(@__DIR__,"heat_transfer_atmosphere.csv"))
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
    linewidth = linewidth,guidefontsize = guidefontsize, dpi = dpi
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

savefig(p7, joinpath(@__DIR__,"heat_transfer_north.png"))

# Diffusion coefficients by surface type

lats = LinRange(-90,90,181)

D_ocean_poles = 0.4
D_ocean_equ = 0.65
D_equ = 0.65
D_NP = 0.28
D_SP = 0.2

function diff_coeff(latitude, type)
    colatitude = deg2rad(90 - latitude) # convert to radians
    if type=="ocean"
        return D_ocean_poles + (D_ocean_equ - D_ocean_poles) * sin(colatitude)^5
    elseif latitude > 0
        return D_NP + (D_equ-D_NP)*sin(colatitude)^5
    else
        return D_SP + (D_equ-D_SP)* sin(colatitude)^5
    end
end

p8 = plot(
    lats,
    diff_coeff.(lats, "ocean"),
    xticks = LinRange(-90,90,7),
    yticks = [0,0.2,0.4,0.6,0.8,1],
    yrange = (0,1),
    xrange = (-90,90), 
    # label = L"\textsf{over} \quad ocean",
    label = "over ocean",
    guidefontsize = guidefontsize,
    xlabel = L"Latitude \; [°]",
    ylabel = L"\widetilde{D} \; [W/m^2/K]",
    size = (700,550),
    linewidth = linewidth,
    dpi = 300
)

plot!(
    p8,
    lats[1:91],
    diff_coeff.(lats[1:91],"not ocean"),
    linewidth = linewidth,
    color = "red",
    label = "over land, snow, sea ice"

)

plot!(
    p8,
    lats[91:end],
    diff_coeff.(lats[91:end],"not ocean"),
    linewidth = linewidth,
    color = "red",
    label =""
)
savefig(p8, joinpath(@__DIR__,"diffusion.png"))