include("milestone1.jl")
include("milestone2.jl")
include("milestone3.jl")
include("milestone4.jl")
include("milestone5.jl")

using SparseArrays

function timestep_euler_backward_2d(jacobian, delta_t)
    A = factorize(sparse(I - delta_t * jacobian))

    function timestep_function(temperature, t, delta_t,
                               mesh, diffusion_coeff, heat_capacity, solar_forcing,
                               radiative_cooling)
        if t == 1
            t_old = size(temperature, 3)
        else
            t_old = t - 1
        end
        # Similar to MS3, we have to solve the equation
        # T_t = T_{t-1} + delta_t * f(T_t, t),
        # where f(T, t) = R(T) + F(t).
        # We use the fact that R is linear and thus can be written as R(T) = AT, where A is the Jacobian of R.
        # Solving for T_t yields
        # T_t = (I - delta_t * A)^{-1} * (T_{t-1} + delta_t * F(t)).
        @views source_terms = calc_source_terms_ebm_2d(heat_capacity,
                                                       solar_forcing[:, :, t],
                                                       radiative_cooling)

        @views temperature[:, :, t] = reshape(A \ vec(temperature[:, :, t_old] +
                                                  delta_t * source_terms),
                                              (mesh.n_latitude, mesh.n_longitude))
    end

    return timestep_function
end

function co2_evolution(jacobian, mesh, diffusion_coeff, heat_capacity, solar_forcing)
    # Read CO2 data
    co2_data = readdlm(joinpath(@__DIR__, "input", "co2_nasa.dat"))

    # Assume that only data for full years is available
    n_years = Int(size(co2_data, 1) / 12)

    average_co2 = [sum(co2_data[(12y + 1):(12y + 12), 4]) / 12 for y in 0:(n_years - 1)]

    ntimesteps = size(solar_forcing, 3)

    average_temperatures = zeros(n_years)
    annual_temperatures = zeros(ntimesteps * n_years)

    timestep_function = timestep_euler_backward_2d(jacobian, 1 / ntimesteps)

    temperature_grid = 15 * ones((mesh.n_latitude, mesh.n_longitude, size(solar_forcing, 3)))

    for y in 1:n_years
        radiative_cooling = calc_radiative_cooling_co2(average_co2[y])
        temperature_grid, area_mean_temp = compute_equilibrium_2d(timestep_function,
                                                                  mesh, diffusion_coeff,
                                                                  heat_capacity,
                                                                  solar_forcing,
                                                                  radiative_cooling,
                                                                  rel_error=1e-2,
                                                                  initial_temperature=temperature_grid)
        annual_temperatures[(ntimesteps * (y - 1) + 1):(ntimesteps * y)] = area_mean_temp
        average_temperatures[y] = sum(area_mean_temp) / ntimesteps
    end

    first_year = Int(co2_data[1, 1])
    last_year = Int(co2_data[end, 1])

    return annual_temperatures, average_temperatures, first_year, last_year
end

function plot_co2_evolution(jacobian, mesh, diffusion_coeff, heat_capacity, solar_forcing)
    annual_temperatures, average_temperatures,
    first_year, last_year = co2_evolution(jacobian, mesh, diffusion_coeff,
                                          heat_capacity, solar_forcing)

    n_timesteps = length(annual_temperatures)

    average_temperatures_per_month = [average_temperatures[floor(Int, (t - 1) / 48) + 1]
                                      for t in 1:n_timesteps]

    labels = first_year:10:last_year

    p = plot(average_temperatures_per_month, label="average temperature",
             xlims=(1, n_timesteps), xticks=(1:480:n_timesteps, labels),
             ylabel="surface temperature [°C]",
             title="Annual temperature with CO2 data from NASA")
    plot!(p, annual_temperatures, label="annual temperature")

    display(p)
end

# Run code
function milestone6()
    geo_dat = read_geography(joinpath(@__DIR__, "input", "The_World128x65.dat"))
    mesh = Mesh(geo_dat)

    albedo = calc_albedo(geo_dat)
    heat_capacity = calc_heat_capacity(geo_dat)

    # Compute solar forcing
    true_longitude = read_true_longitude(joinpath(@__DIR__, "input", "True_Longitude.dat"))
    solar_forcing = calc_solar_forcing(albedo, true_longitude)
    ntimesteps = length(true_longitude)

    # Compute and plot diffusion coefficient
    diffusion_coeff = calc_diffusion_coefficients(geo_dat)

    co2_ppm = 315.0
    radiative_cooling = calc_radiative_cooling_co2(co2_ppm)

    jacobian = calc_jacobian_ebm_2d(mesh, diffusion_coeff, heat_capacity)

    temperature, _ = compute_equilibrium_2d(timestep_euler_backward_2d(jacobian,
                                                                       1 / ntimesteps),
                                            mesh, diffusion_coeff, heat_capacity,
                                            solar_forcing,
                                            radiative_cooling)

    plot_temperature(temperature, geo_dat, 1)

    # Copied from MS4
    # Area mean of pointwise annual temperature
    annual_mean_temperature_north = [calc_mean_north(temperature[:, :, t], mesh.area)
                                     for t in 1:ntimesteps]
    annual_mean_temperature_south = [calc_mean_south(temperature[:, :, t], mesh.area)
                                     for t in 1:ntimesteps]
    annual_mean_temperature_total = [calc_mean(temperature[:, :, t], mesh.area)
                                     for t in 1:ntimesteps]

    average_temperature_north = sum(annual_mean_temperature_north) / ntimesteps
    average_temperature_south = sum(annual_mean_temperature_south) / ntimesteps
    average_temperature_total = sum(annual_mean_temperature_total) / ntimesteps

    plot_annual_temperature_north_south(annual_mean_temperature_north,
                                        annual_mean_temperature_south,
                                        annual_mean_temperature_total,
                                        average_temperature_north,
                                        average_temperature_south,
                                        average_temperature_total)

    # Compute temperature in Cologne.
    # Cologne lies about halfway between these two grid points.
    annual_temperature_cologne = (temperature[15, 68, :] +
                                  temperature[15, 69, :]) / 2
    average_temperature_cologne = sum(annual_temperature_cologne) / ntimesteps

    plot_annual_temperature(annual_temperature_cologne, average_temperature_cologne,
                            "Annual temperature with CO2 = $co2_ppm [ppm] in Cologne")

    plot_co2_evolution(jacobian, mesh, diffusion_coeff, heat_capacity, solar_forcing)

    # Animate annual temperature
    anim = @animate for ts in 1:ntimesteps
        plot_temperature(temperature, geo_dat, ts)
    end

    gif(anim, joinpath(@__DIR__, "annual_temperature.gif"), fps=7)
end
