using PlotlyJS
using DelimitedFiles
using LaTeXStrings

include("milestone6.jl")


function parametrisation(surface_data, radius=10)

    lat_resolution, long_resolution = size(surface_data)

    latitude = range(0, pi, lat_resolution)
    longitude = range(0, 2*pi, long_resolution)

    lat_grid = latitude'.*ones(long_resolution)
    long_grid = ones(lat_resolution)'.*longitude

    X = radius*sin.(lat_grid).*cos.(long_grid)
    Y = radius*sin.(lat_grid).*sin.(long_grid)
    Z = radius*cos.(lat_grid)

    return X,Y,Z

end

function plot_earth(X,Y,Z,surface_data)

    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = surface_data',
        colorscale = [
            [0, "rgb(60,149,33)"], #land, green
            [0.25, "rgb(165,255,255)"], #sea ice, icy blue
            [0.5, "rgb(255,255,255)"], #snow, white
            [1, "rgb(26,100,212)" ], #ocean, blue
        ],
        colorbar = (
            autotick = false, 
            tick0 = 0,
            dtick = 1,
            tickcolor = 888,
            tickvals = [1,2,3,4,5],
            ticktext = [
                "land",
                "sea ice",
                "snow cover",
                "lakes",
                "ocean"
            ],
            tickfont = (
                color = "rgb(255,255,255)",
                size = 15
            )
        )
    )
    
    layout = Layout(
                scene = attr(
                    xaxis = attr(
                        visible = false
                    ),
                    yaxis = attr(
                        visible = false
                    ),
                    zaxis = attr(
                        visible = false
                    ),
                ),
                paper_bgcolor = "black",
                title = attr(
                    text = "The Earth",
                    y=0.95,
                    x=0.5
                ),
                titlefont = (
                    color = "rgb(255,255,255)",
                    size = 40
                ),
    )

    return Plot(fig1, layout)   
end

function get_outlines(geo_dat)
    outlines = zeros(size(geo_dat))
    for i in 2:64 # Edge cases can be ignored because there is no continent 
        for j in 2:127
            if geo[i,j] == 1
                if ((geo[i,j-1] != 1) || (geo[i,j+1] != 1) || (geo[i-1,j] != 1) || (geo[i+1,j] != 1))
                    outlines[i,j] = 1
                end
            end
        end
    end
    return round.(Int, outlines)
end

function plot_albedo_3d(X,Y,Z, albedo, surface_data)

    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = surface_data',
        showscale = false,
        colorscale = [
            [0, "rgb(255,255,255)"], # no border
            [1,"rgb(0,0,0)" ], # continent border
        ],
    )
    
    fig2 = PlotlyJS.surface(
        x = 1.01 .* X,
        y = 1.01 .* Y,
        z = 1.01 .* Z,
        surfacecolor = albedo',
        showscale = true,
        opacity = 0.5,
        colorscale = [
            [0, "rgb(0,0,0)"], # black, nothing reflected
            [1,"rgb(255,255,255)" ], # white, everything reflected
        ],
        colorbar = (
            autotick = false, 
            tickcolor = 888,
            tickfont = (
                color = "rgb(255,255,255)",
                size = 15
                ),
                
            title="albedo",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 15
            )  
        )
    )
    
    layout = Layout(
                scene = attr(
                    xaxis = attr(
                        visible = false
                    ),
                    yaxis = attr(
                        visible = false
                    ),
                    zaxis = attr(
                        visible = false
                    ),
                ),
                paper_bgcolor = "black",
                title = attr(
                    text = "Surface albedo",
                    y=0.95,
                    x=0.5
                ),
                titlefont = (
                    color = "rgb(255,255,255)",
                    size = 40
                ),
            )


    return Plot([fig1,fig2], layout)   
end

function plot_heatcapacity_3d(X,Y,Z, heat_capacity, surface_data)

    heat_capacity_log = log10.(heat_capacity)

    vmin = minimum(heat_capacity_log)
    vmax = maximum(heat_capacity_log)

    cbar_ticks = LinRange(vmin, vmax, 10)
    cbar_tick_values = [@sprintf("%.3f", 10^tick) for tick in cbar_ticks]

    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = surface_data',
        showscale = false,
        colorscale = [
            [0, "rgb(255,255,255)"], # no border
            [1,"rgb(0,0,0)" ], # continent border
        ],
    )
    
    fig2 = PlotlyJS.surface(
        x = 1.01 .* X,
        y = 1.01 .* Y,
        z = 1.01 .* Z,
        surfacecolor = heat_capacity_log',
        showscale = true,
        opacity = 0.7,
        colorscale =[
            [0, "rgb(255,255,255)" ], # white
            [1, "rgb(139,0,0)"], # red 
        ],
        colorbar = (
            autotick = false, 
            tickcolor = 888,
            tickfont = (
                color = "rgb(255,255,255)",
                size = 15
                ),
            tickvals = cbar_ticks,
            ticktext = cbar_tick_values,
            title="[ J / m² / K ]",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 20
                ),
            titleside = "right"   
        )
    )
    
    layout = Layout(
                scene = attr(
                    xaxis = attr(
                        visible = false
                    ),
                    yaxis = attr(
                        visible = false
                    ),
                    zaxis = attr(
                        visible = false
                    ),
                ),
                paper_bgcolor = "black",
                title = attr(
                    text = "Heat capacity",
                    y=0.95,
                    x=0.5
                ),
                titlefont = (
                    color = "rgb(255,255,255)",
                    size = 40
                ),
    )

    return Plot([fig1,fig2], layout)   
end

function plot_solar_forcing_3d_anim(X,Y,Z,solar_forcing, surface_data)

    cmin, cmax = extrema(solar_forcing)
    
    n_frames = size(solar_forcing, 3)
    
    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = get_outlines(surface_data)',
        showscale = false,
        colorscale = [
            [0, "rgb(255,255,255)"], # no border
            [1,"rgb(0,0,0)" ], # continent border
        ],
    )

    fig2 = PlotlyJS.surface(
        x = 1.01 .* X,
        y = 1.01 .* Y,
        z = 1.01 .* Z,
        surfacecolor = solar_forcing[:,:,1]',
        showscale = true,
        opacity = 0.8,
        colorscale = "Hot",
        cmax = cmax,
        cmin = cmin,
        colorbar = (
            autotick = false, 
            tickcolor = 888,
            tickfont = (
                color = "rgb(255,255,255)",
                size = 20
            ),
            title="[ W / m² ]",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 20
            ),
            titleside = "right"   
        )
    )

    

    sliders = [
        attr(
            steps = [
                attr(
                    method = "animate",
                    args= [
                        [
                            "fr$k"
                        ],                           
                        attr(
                            mode = "immediate",
                            frame = attr(
                                duration = 40,
                                redraw = true
                            ),
                            transition = attr(
                                duration = 0
                            )
                        )
                    ],
                    label = "$k"
                )
            for k in 1:n_frames], 
            active = 17,
            transition = attr(
                duration = 0
            ),
            x=0, # slider starting position  
            y=0, 
            currentvalue = attr(
                font = attr(
                        size = 12
                    ), 
                prefix = "Step: ", 
                visible = true, 
                xanchor = "center"
            ),  
            len = 1.0 #slider length
        )    
    ];
    
    updatemenus = [
        attr(
            type = "buttons", 
            active = 0,
            y = 0.0,  #(x,y) button position 
            x = 1,
            buttons = [
                attr(
                    label = "Play",
                    method = "animate",
                    args = [
                        nothing,
                        attr(
                            frame = attr(
                                duration = 5, 
                                redraw = true
                            ),
                            transition = attr(
                                duration = 0
                            ),
                            fromcurrent = true,
                            mode = "immediate"
                        )
                    ]
                )
            ]
        )
    ];



    frames  = Vector{PlotlyFrame}(undef, n_frames)
    for k in 1:n_frames
        day = (round(Int, (k - 1) / n_frames * 365) + 80) % 365
        frames[k] = PlotlyJS.frame(
                        data = [
                            attr(
                                surfacecolor = solar_forcing[:,:,k]',            
                                colorbar = attr(
                                    #title = "Day $day",
                                )
                            )
                        ],
                        layout = attr(
                            title_text = "Solar forcing - Day $day"
                        ),
                        name="fr$k",
                        traces=[1],
                    )
    end    

    layout = Layout(
        scene = attr(
            xaxis = attr(
                visible = false
            ),
            yaxis = attr(
                visible = false
            ),
            zaxis = attr(
                visible = false
            ),
        ),
        paper_bgcolor = "black",
        sliders = sliders,
        title = attr(
            text = "Solar forcing",
            y=0.95,
            x=0.5
        ),
        titlefont = (
            color = "rgb(255,255,255)",
            size = 40
            ),
        updatemenus = updatemenus
    )

    return Plot([fig1,fig2], layout, frames)
end

function plot_diffusioncoefficient_3d(X,Y,Z, diffusion_coefficient, surface_data)

    # vmin = minimum(diffusion_coefficient)
    # vmax = maximum(diffusion_coefficient)

    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = surface_data',
        showscale = false,
        colorscale = [
            [0, "rgb(255,255,255)"], # no border
            [1,"rgb(0,0,0)" ], # continent border
        ],
    )
    
    fig2 = PlotlyJS.surface(
        x = 1.01 .* X,
        y = 1.01 .* Y,
        z = 1.01 .* Z,
        surfacecolor = diffusion_coefficient',
        showscale = true,
        opacity = 0.75,
        colorscale =[
            [0, "rgb(0,0,255)" ], # blue
            [1, "rgb(255,255,0)"], # yellow 
        ],
        colorbar = (
            autotick = false, 
            tickcolor = 888,
            tickfont = (
                color = "rgb(255,255,255)",
                size = 15
                ),
            title="[ W / K ]",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 20
                ),
            titleside = "right"   
        )
    )
    
    layout = Layout(
                scene = attr(
                    xaxis = attr(
                        visible = false
                    ),
                    yaxis = attr(
                        visible = false
                    ),
                    zaxis = attr(
                        visible = false
                    ),
                ),
                paper_bgcolor = "black",
                title = attr(
                    text = "Diffusion coefficient of the 2D EBM",
                    x = 0.5,
                    y = 0.95
                ),
                titlefont = (
                    color = "rgb(255,255,255)",
                    size = 40
                    ),

            )


    return Plot([fig1,fig2], layout)   
end

function annual_temperature_pointwise(nlatitude, nlongitude, ntimesteps, heat_capacity, solar_forcing, radiative_cooling)
    # Calculate annual temperature for every grid point
    annual_temperature_pointwise = Array{Float64, 3}(undef, nlatitude, nlongitude,
                                                     ntimesteps)
    for j in 1:nlongitude, i in 1:nlatitude
        # I/O in Julia is pretty slow, so this takes forever with `verbose=true`.
        annual_temperature, _ = compute_equilibrium(timestep_euler_forward,
                                                    heat_capacity[i, j],
                                                    solar_forcing[i, j, :],
                                                    radiative_cooling,
                                                    verbose=false)

        annual_temperature_pointwise[i, j, :] = annual_temperature
    end
    return annual_temperature_pointwise
end

function plot_temperature_3d_anim(X,Y,Z,temp, surface_data)

    cmin, cmax = extrema(temp)
    if (abs(cmax)> abs(cmin))
        cmin = -cmax
    else
        cmax = -cmin
    end

    n_frames = size(temp, 3)
    
    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = get_outlines(surface_data)',
        showscale = false,
        colorscale = [
            [0, "rgb(255,255,255)"], # no border
            [1,"rgb(0,0,0)" ], # continent border
        ],
    )

    fig2 = PlotlyJS.surface(
        x = 1.01 .* X,
        y = 1.01 .* Y,
        z = 1.01 .* Z,
        surfacecolor = temp[:,:,1]',
        showscale = true,
        opacity = 0.8,
        colorscale = "RdBu",
        cmax = cmax,
        cmin = cmin,
        colorbar = (
            autotick = false, 
            tickcolor = 888,
            tickfont = (
                color = "rgb(255,255,255)",
                size = 20
            ),
            title="[ ° C ]",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 20
            ),
            titleside = "right"   
        )
    )

    

    sliders = [
        attr(
            steps = [
                attr(
                    method = "animate",
                    args= [
                        [
                            "fr$k"
                        ],                           
                        attr(
                            mode = "immediate",
                            frame = attr(
                                duration = 40,
                                redraw = true
                            ),
                            transition = attr(
                                duration = 0
                            )
                        )
                    ],
                    label = "$k"
                )
            for k in 1:n_frames], 
            active = 17,
            transition = attr(
                duration = 0
            ),
            x=0, # slider starting position  
            y=0, 
            currentvalue = attr(
                font = attr(
                        size = 12
                    ), 
                prefix = "Step: ", 
                visible = true, 
                xanchor = "center"
            ),  
            len = 1.0 #slider length
        )    
    ];
    
    updatemenus = [
        attr(
            type = "buttons", 
            active = 0,
            y = 0.0,  #(x,y) button position 
            x = 1,
            buttons = [
                attr(
                    label = "Play",
                    method = "animate",
                    args = [
                        nothing,
                        attr(
                            frame = attr(
                                duration = 5, 
                                redraw = true
                            ),
                            transition = attr(
                                duration = 0
                            ),
                            fromcurrent = true,
                            mode = "immediate"
                        )
                    ]
                )
            ]
        )
    ];



    frames  = Vector{PlotlyFrame}(undef, n_frames)
    for k in 1:n_frames
        day = (round(Int, (k - 1) / n_frames * 365) + 80) % 365
        frames[k] = PlotlyJS.frame(
                        data = [
                            attr(
                                surfacecolor = temp[:,:,k]',            
                                colorbar = attr(
                                    #title = "Day $day",
                                )
                            )
                        ],
                        layout = attr(
                            title_text = "Temperature - Day $day"
                        ),
                        name="fr$k",
                        traces=[1],
                    )
    end    

    layout = Layout(
        scene = attr(
            xaxis = attr(
                visible = false
            ),
            yaxis = attr(
                visible = false
            ),
            zaxis = attr(
                visible = false
            ),
        ),
        paper_bgcolor = "black",
        sliders = sliders,
        title = attr(
            text = "Temperature",
            y=0.95,
            x=0.5
        ),
        titlefont = (
            color = "rgb(255,255,255)",
            size = 40
            ),
        updatemenus = updatemenus
    )

    return Plot([fig1,fig2], layout, frames)
end



geo = readdlm("The_World128x65.dat")

albedo = calc_albedo(geo)
heat_capacity = calc_heat_capacity(geo)
true_lon = read_true_longitude("True_Longitude.dat")
solar_forcing = calc_solar_forcing(albedo,true_lon)
X,Y,Z = parametrisation(geo)
area = calc_area(geo) # Compute area-mean quantities
outlines = get_outlines(geo)
co2_ppm = 315.0
radiative_cooling = calc_radiative_cooling_co2(co2_ppm)

nlatitude, nlongitude = size(geo)
ntimesteps = length(true_lon)

#temperature = annual_temperature_pointwise(nlatitude, nlongitude,ntimesteps,heat_capacity,solar_forcing,radiative_cooling)
diffusion_coefficient = calc_diffusion_coefficients(geo)
diffusion_coeff = diffusion_coefficient

mesh = Mesh(geo)

jacobian = calc_jacobian_ebm_2d(mesh, diffusion_coeff, heat_capacity)

temperature,_ = compute_equilibrium_2d(timestep_euler_backward_2d(jacobian, 1/48), mesh, diffusion_coefficient, heat_capacity, solar_forcing, radiative_cooling)


# display(plot_earth(X,Y,Z,geo)) # earth
# display(plot_albedo_3d(X,Y,Z, albedo, outlines))# albedo
# display(plot_heatcapacity_3d(X,Y,Z, heat_capacity, outlines)) # heat capacity
# display(plot_diffusioncoefficient_3d(X,Y,Z, diffusion_coefficient, outlines))


# display(plot_solar_forcing_3d_anim(X,Y,Z,solar_forcing,outlines))

# display(plot_temperature_3d_anim(X,Y,Z,temperature,outlines))

earth = plot_earth(X,Y,Z,geo)
open("./earth.html", "w") do io
    PlotlyBase.to_html(io, earth)
end

albedo_plot = plot_albedo_3d(X,Y,Z, albedo, outlines)
# Plots.savefig(albedo_plot, "savealb.html")
open("./albedo.html", "w") do io
    PlotlyBase.to_html(io, albedo_plot)
end

solar_forcing_plot = plot_solar_forcing_3d_anim(X,Y,Z,solar_forcing, outlines)
open("./sf.html", "w") do io
    PlotlyBase.to_html(io, solar_forcing_plot)
end

heat_capacity_plot = plot_heatcapacity_3d(X,Y,Z, heat_capacity, outlines)
open("./heat_capacity.html", "w") do io
    PlotlyBase.to_html(io, heat_capacity_plot)
end

diff_coeff_plot = plot_diffusioncoefficient_3d(X,Y,Z, diffusion_coefficient, outlines)
open("./diffusion_coeff.html", "w") do io
    PlotlyBase.to_html(io, diff_coeff_plot)
end

temperature_plot = plot_temperature_3d_anim(X,Y,Z,temperature,outlines)
open("./temperature.html", "w") do io
    PlotlyBase.to_html(io, temperature_plot)
end