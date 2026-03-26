include("milestone1.jl")
include("milestone2.jl")
include("milestone3.jl")

function calc_mean_north(data, area)
    nlatitude, nlongitude = size(data)
    j_equator = Int((nlatitude - 1) / 2) + 1

    # North Pole
    mean_data = area[1] * data[1, 1]

    # Inner nodes
    for j in 2:(j_equator - 1)
        for i in 1:nlongitude
            mean_data += area[j] * data[j, i]
        end
    end

    # Equator
    for i in 1:nlongitude
        mean_data += 0.5 * area[j_equator] * data[j_equator, i]
    end

    return 2 * mean_data
end

function calc_mean_south(data, area)
    nlatitude, nlongitude = size(data)
    j_equator = Int((nlatitude - 1) / 2) + 1

    # South Pole
    mean_data = area[end] * data[end, end]

    # Inner nodes
    for j in (j_equator + 1):(nlatitude - 1)
        for i in 1:nlongitude
            mean_data += area[j] * data[j, i]
        end
    end

    # Equator
    for i in 1:nlongitude
        mean_data += 0.5 * area[j_equator] * data[j_equator, i]
    end

    return 2 * mean_data
end

function plot_annual_temperature_north_south(annual_temperature_north,
                                             annual_temperature_south,
                                             annual_temperature_total,
                                             average_temperature_north,
                                             average_temperature_south,
                                             average_temperature_total)
    ntimesteps = length(annual_temperature_total)
    labels = ["March", "June", "September", "December", "March"]

    p = plot(average_temperature_total * ones(ntimesteps),
             label="average temperature (total)",
             xlims=(1, ntimesteps), xticks=(LinRange(1, ntimesteps, 5), labels),
             ylabel="surface temperature [°C]",
             title="Annual temperature with CO2 = 315 [ppm]")
    plot!(p, average_temperature_north * ones(ntimesteps),
          label="average temperature (north)")
    plot!(p, average_temperature_south * ones(ntimesteps),
          label="average temperature (south)")
    plot!(p, annual_temperature_total, label="temperature (total)")
    plot!(p, annual_temperature_north, label="temperature (north)")
    plot!(p, annual_temperature_south, label="temperature (south)")

    display(p)
end

function plot_temperature(temperature, geo_dat, timestep)
    vmin = minimum(temperature)
    vmax = -vmin # We want to have 0°C in the center.

    nlatitude, nlongitude = size(temperature)
    x, y = robinson_projection(nlatitude, nlongitude)

    ntimesteps = size(temperature, 3)
    day = (round(Int, (timestep - 1) / ntimesteps * 365) + 80) % 365
    plot = contourf(x, y, temperature[:, :, timestep],
                    clims=(vmin, vmax),
                    levels=LinRange(vmin, vmax, 200),
                    aspect_ratio=1,
                    title="Temperature for Day $day",
                    c=:seismic,
                    colorbar_title="temperature [°C]",
                    axis=([], false),
                    dpi=300)

    # Add contour lines
    contour!(plot, x, y, geo_dat, levels=[0.5, 1.5, 2.5, 3.5, 6.5],
             color=[:black], linewidth=0.6)

    return plot
end

# Run code
function milestone4()
    geo_dat = read_geography(joinpath(@__DIR__, "input", "The_World128x65.dat"))
    nlatitude, nlongitude = size(geo_dat)

    albedo = calc_albedo(geo_dat)
    heat_capacity = calc_heat_capacity(geo_dat)

    # Compute solar forcing
    true_longitude = read_true_longitude(joinpath(@__DIR__, "input", "True_Longitude.dat"))
    solar_forcing = calc_solar_forcing(albedo, true_longitude)

    # Compute area-mean quantities
    area = calc_area(geo_dat)

    mean_albedo_north = calc_mean_north(albedo, area)
    print("Mean albedo north = $mean_albedo_north")
    mean_albedo_south = calc_mean_south(albedo, area)
    print("Mean albedo south = $mean_albedo_south")
    mean_albedo_total = calc_mean(albedo, area)
    print("Mean albedo total = $mean_albedo_total")

    mean_heat_capacity_north = calc_mean_north(heat_capacity, area)
    print("Mean heat capacity north = $mean_heat_capacity_north")
    mean_heat_capacity_south = calc_mean_south(heat_capacity, area)
    print("Mean heat capacity south = $mean_heat_capacity_south")
    mean_heat_capacity_total = calc_mean(heat_capacity, area)
    print("Mean heat capacity total = $mean_heat_capacity_total")

    ntimesteps = length(true_longitude)
    mean_solar_forcing_north = [calc_mean_north(solar_forcing[:, :, t], area)
                                for t in 1:ntimesteps]
    mean_solar_forcing_south = [calc_mean_south(solar_forcing[:, :, t], area)
                                for t in 1:ntimesteps]
    mean_solar_forcing_total = [calc_mean(solar_forcing[:, :, t], area)
                                for t in 1:ntimesteps]

    co2_ppm = 315.0
    radiative_cooling = calc_radiative_cooling_co2(co2_ppm)

    # Compute equilibrium for all three means
    annual_temperature_north_, average_temperature_north_ = compute_equilibrium(timestep_euler_forward,
                                                                                mean_heat_capacity_north,
                                                                                mean_solar_forcing_north,
                                                                                radiative_cooling)
    annual_temperature_south_, average_temperature_south_ = compute_equilibrium(timestep_euler_forward,
                                                                                mean_heat_capacity_south,
                                                                                mean_solar_forcing_south,
                                                                                radiative_cooling)
    annual_temperature_total_, average_temperature_total_ = compute_equilibrium(timestep_euler_forward,
                                                                                mean_heat_capacity_total,
                                                                                mean_solar_forcing_total,
                                                                                radiative_cooling)

    plot_annual_temperature_north_south(annual_temperature_north_,
                                        annual_temperature_south_,
                                        annual_temperature_total_,
                                        average_temperature_north_,
                                        average_temperature_south_,
                                        average_temperature_total_)

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

    # Area mean of pointwise annual temperature
    annual_mean_temperature_north = [calc_mean_north(annual_temperature_pointwise[:, :, t],
                                                     area)
                                     for t in 1:ntimesteps]
    annual_mean_temperature_south = [calc_mean_south(annual_temperature_pointwise[:, :, t],
                                                     area)
                                     for t in 1:ntimesteps]
    annual_mean_temperature_total = [calc_mean(annual_temperature_pointwise[:, :, t], area)
                                     for t in 1:ntimesteps]

    average_temperature_north_ = sum(annual_mean_temperature_north) / ntimesteps
    average_temperature_south_ = sum(annual_mean_temperature_south) / ntimesteps
    average_temperature_total_ = sum(annual_mean_temperature_total) / ntimesteps

    plot_annual_temperature_north_south(annual_mean_temperature_north,
                                        annual_mean_temperature_south,
                                        annual_mean_temperature_total,
                                        average_temperature_north_,
                                        average_temperature_south_,
                                        average_temperature_total_)

    # Compute temperature in Cologne.
    # Cologne lies about halfway between these two grid points.
    annual_temperature_cologne = (annual_temperature_pointwise[15, 68, :] +
                                  annual_temperature_pointwise[15, 69, :]) / 2
    average_temperature_cologne = sum(annual_temperature_cologne) / ntimesteps

    plot_annual_temperature(annual_temperature_cologne, average_temperature_cologne,
                            "Annual temperature with CO2 = $co2_ppm [ppm] in Cologne")

    # Animate annual temperature
    anim = @animate for ts in 1:ntimesteps
        plot_temperature(annual_temperature_pointwise, geo_dat, ts)
    end

    gif(anim, joinpath(@__DIR__, "annual_temperature.gif"), fps=7)
end
