import os

import imageio
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

from milestone1 import read_geography, robinson_projection
from milestone2 import read_true_longitude
from milestone3 import calc_area, calc_albedo, calc_solar_forcing, calc_mean, calc_heat_capacity, \
    calc_radiative_cooling_co2, compute_equilibrium, timestep_euler_forward, plot_annual_temperature


def calc_mean_north(data, area):
    nlatitude, nlongitude = data.shape
    j_equator = int((nlatitude - 1) / 2)

    # North Pole
    mean_data = area[0] * data[0, 0]

    # Inner nodes
    for j in range(1, j_equator):
        for i in range(nlongitude):
            mean_data += area[j] * data[j, i]

    # Equator
    for i in range(nlongitude):
        mean_data += 0.5 * area[j_equator] * data[j_equator, i]

    return 2 * mean_data


def calc_mean_south(data, area):
    nlatitude, nlongitude = data.shape
    j_equator = int((nlatitude - 1) / 2)

    # South Pole
    mean_data = area[-1] * data[-1, -1]

    # Inner nodes
    for j in range(j_equator + 1, nlatitude - 1):
        for i in range(nlongitude):
            mean_data += area[j] * data[j, i]

    # Equator
    for i in range(nlongitude):
        mean_data += 0.5 * area[j_equator] * data[j_equator, i]

    return 2 * mean_data


def plot_annual_temperature_north_south(annual_temperature_north, annual_temperature_south, annual_temperature_total,
                                        average_temperature_north, average_temperature_south,
                                        average_temperature_total):
    fig, ax = plt.subplots()

    ntimesteps = len(annual_temperature_total)
    plt.plot(average_temperature_total * np.ones(ntimesteps), label="average temperature (total)")
    plt.plot(average_temperature_north * np.ones(ntimesteps), label="average temperature (north)")
    plt.plot(average_temperature_south * np.ones(ntimesteps), label="average temperature (south)")
    plt.plot(annual_temperature_total, label="temperature (total)")
    plt.plot(annual_temperature_north, label="temperature (north)")
    plt.plot(annual_temperature_south, label="temperature (south)")

    plt.xlim((0, ntimesteps - 1))
    labels = ["March", "June", "September", "December", "March"]
    plt.xticks(np.linspace(0, ntimesteps - 1, 5), labels)
    ax.set_ylabel("surface temperature [°C]")
    plt.grid()
    plt.title(f"Annual temperature with CO2 = 315 [ppm]")
    plt.legend(loc="upper right")

    plt.tight_layout()
    plt.show()


# Similar to plot_robinson_projection in MS1, but we also add contour lines
def plot_robinson_projection_with_lines(data, geo_dat, title, **kwargs):
    # Get the coordinates for the Robinson projection.
    nlatitude, nlongitude = data.shape
    x, y = robinson_projection(nlatitude, nlongitude)

    # Start plotting.
    fig, ax = plt.subplots()

    # Create contour plot of geography information against x and y.
    im = ax.contourf(x, y, data, **kwargs)

    # Add contour lines with levels from MS1
    levels = [0.5, 1.5, 2.5, 3.5, 6.5]
    ax.contour(x, y, geo_dat, colors="black", linewidths=0.6, levels=levels, linestyles="solid")

    plt.title(title)
    ax.set_aspect("equal")

    # Remove axes and ticks.
    plt.xticks([])
    plt.yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Colorbar with the same height as the plot. Code copied from
    # https://stackoverflow.com/a/18195921
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)

    return cbar


def plot_temperature(temperature, geo_dat, timestep, show_plot=False):
    vmin = np.amin(temperature)
    vmax = np.amax(temperature)
    levels = np.linspace(vmin, vmax, 200)
    norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    # Reuse plotting function from milestone 1.
    ntimesteps = temperature.shape[2]
    day = (np.int_(timestep / ntimesteps * 365) + 80) % 365
    cbar = plot_robinson_projection_with_lines(temperature[:, :, timestep], geo_dat,
                                               f"Temperature for Day {day}",
                                               levels=levels, cmap="seismic",
                                               vmin=vmin, vmax=vmax, norm=norm)
    cbar.set_label("surface temperature [°C]")

    # Adjust size of plot to viewport to prevent clipping of the legend.
    plt.tight_layout()

    filename = f"temperature_{timestep}.png"
    plt.savefig(filename, dpi=300)
    if show_plot:
        plt.show()
    plt.close()

    return filename


# Run code
if __name__ == '__main__':
    geo_dat_ = read_geography("input/The_World128x65.dat")
    nlatitude_, nlongitude_ = geo_dat_.shape

    albedo = calc_albedo(geo_dat_)
    heat_capacity = calc_heat_capacity(geo_dat_)

    # Compute solar forcing
    true_longitude = read_true_longitude("input/True_Longitude.dat")
    solar_forcing = calc_solar_forcing(albedo, true_longitude)

    # Compute area-mean quantities
    area_ = calc_area(geo_dat_)

    mean_albedo_north = calc_mean_north(albedo, area_)
    print(f"Mean albedo north = {mean_albedo_north}")
    mean_albedo_south = calc_mean_south(albedo, area_)
    print(f"Mean albedo south = {mean_albedo_south}")
    mean_albedo_total = calc_mean(albedo, area_)
    print(f"Mean albedo total = {mean_albedo_total}")

    mean_heat_capacity_north = calc_mean_north(heat_capacity, area_)
    print(f"Mean heat capacity north = {mean_heat_capacity_north}")
    mean_heat_capacity_south = calc_mean_south(heat_capacity, area_)
    print(f"Mean heat capacity south = {mean_heat_capacity_south}")
    mean_heat_capacity_total = calc_mean(heat_capacity, area_)
    print(f"Mean heat capacity total = {mean_heat_capacity_total}")

    ntimesteps_ = len(true_longitude)
    mean_solar_forcing_north = [calc_mean_north(solar_forcing[:, :, t], area_) for t in range(ntimesteps_)]
    mean_solar_forcing_south = [calc_mean_south(solar_forcing[:, :, t], area_) for t in range(ntimesteps_)]
    mean_solar_forcing_total = [calc_mean(solar_forcing[:, :, t], area_) for t in range(ntimesteps_)]

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

    plot_annual_temperature_north_south(annual_temperature_north_, annual_temperature_south_, annual_temperature_total_,
                                        average_temperature_north_, average_temperature_south_,
                                        average_temperature_total_)

    # Calculate annual temperature for every grid point
    def annual_temp(i, j):
        annual_temperature, _ = compute_equilibrium(timestep_euler_forward, heat_capacity[i, j],
                                                    solar_forcing[i, j, :], radiative_cooling,
                                                    verbose=False)

        return annual_temperature

    annual_temperature_pointwise = np.array(
        [[annual_temp(i, j) for j in range(nlongitude_)] for i in range(nlatitude_)])

    # Area mean of pointwise annual temperature
    annual_mean_temperature_north = [calc_mean_north(annual_temperature_pointwise[:, :, t], area_)
                                     for t in range(ntimesteps_)]
    annual_mean_temperature_south = [calc_mean_south(annual_temperature_pointwise[:, :, t], area_)
                                     for t in range(ntimesteps_)]
    annual_mean_temperature_total = [calc_mean(annual_temperature_pointwise[:, :, t], area_)
                                     for t in range(ntimesteps_)]

    average_temperature_north_ = np.sum(annual_mean_temperature_north) / ntimesteps_
    average_temperature_south_ = np.sum(annual_mean_temperature_south) / ntimesteps_
    average_temperature_total_ = np.sum(annual_mean_temperature_total) / ntimesteps_

    plot_annual_temperature_north_south(annual_mean_temperature_north, annual_mean_temperature_south,
                                        annual_mean_temperature_total, average_temperature_north_,
                                        average_temperature_south_, average_temperature_total_)

    # Compute temperature in Cologne.
    # Cologne lies about halfway between these two grid points.
    annual_temperature_cologne = (annual_temperature_pointwise[14, 67, :] +
                                  annual_temperature_pointwise[14, 68, :]) / 2
    average_temperature_cologne = np.sum(annual_temperature_cologne) / ntimesteps_

    plot_annual_temperature(annual_temperature_cologne, average_temperature_cologne,
                            f"Annual temperature with CO2 = {co2_ppm} [ppm] in Cologne")

    # Plot temperature for each time step
    filenames = []
    for ts in range(ntimesteps_):
        filename_ = plot_temperature(annual_temperature_pointwise, geo_dat_, ts, show_plot=False)
        filenames.append(filename_)

    # Build GIF
    frames = [imageio.v3.imread(filename_) for filename_ in filenames]
    imageio.mimsave("annual_temperature.gif", frames)

    # Remove files
    for filename_ in set(filenames):
        os.remove(filename_)
