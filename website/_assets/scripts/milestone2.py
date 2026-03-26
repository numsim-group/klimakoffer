import os

import imageio
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter

from milestone1 import read_geography, plot_robinson_projection


def calc_albedo(geo_dat):
    def legendre(latitude):
        return 0.5 * (3 * np.sin(latitude) ** 2 - 1)

    def albedo(surface_type, latitude):
        if surface_type == 1:
            return 0.3 + 0.12 * legendre(latitude)
        elif surface_type == 2:
            return 0.6
        elif surface_type == 3:
            return 0.75
        elif surface_type == 5:
            return 0.29 + 0.12 * legendre(latitude)
        else:
            raise ValueError(f"Unknown surface type {surface_type}.")

    nlatitude, nlongitude = geo_dat.shape
    y_lat = np.linspace(np.pi / 2, -np.pi / 2, nlatitude)

    # Map surface type to albedo.
    return np.array([[albedo(geo_dat[i, j], y_lat[i])
                      for j in range(nlongitude)]
                     for i in range(nlatitude)])


def plot_albedo(albedo):
    # Minimum and maximum of the values of the albedo.
    vmin = 0.05  # np.amin(albedo)
    vmax = np.amax(albedo)

    levels = np.linspace(vmin, vmax, 100)

    # Reuse plotting function from milestone 1.
    cbar = plot_robinson_projection(albedo, "Surface Albedo of the Earth",
                                    levels=levels, cmap="gray",
                                    vmin=vmin, vmax=vmax)
    cbar.set_ticks([vmin, 0.5 * (vmin + vmax), vmax])
    cbar.set_label("albedo")

    # Adjust size of plot to viewport to prevent clipping of the legend.
    plt.tight_layout()
    plt.show()


def calc_heat_capacity(geo_dat):
    sec_per_yr = 3.15576e7  # seconds per year

    c_atmos = 1.225 * 1000 * 3850
    c_ocean = 1030 * 4000 * 70
    c_seaice = 917 * 2000 * 1.5
    c_land = 1350 * 750 * 1
    c_snow = 400 * 880 * 0.5

    def heat_capacity(surface_type):
        if surface_type == 1:
            capacity_surface = c_land
        elif surface_type == 2:
            capacity_surface = c_seaice
        elif surface_type == 3:
            capacity_surface = c_snow
        elif surface_type == 5:
            capacity_surface = c_ocean
        else:
            raise ValueError(f"Unknown surface type {surface_type}.")

        return (capacity_surface + c_atmos) / sec_per_yr

    # Map surface type to heat capacity.
    nlatitude, nlongitude = geo_dat.shape
    return np.array([[heat_capacity(geo_dat[i, j])
                      for j in range(nlongitude)]
                     for i in range(nlatitude)])


def plot_heat_capacity(heat_capacity):
    vmin = np.amin(heat_capacity)
    vmax = np.amax(heat_capacity)

    levels = np.power(10, np.linspace(np.log10(vmin), np.log10(vmax), 500))

    # Reuse plotting function from milestone 1.
    cbar = plot_robinson_projection(heat_capacity, "Surface Heat Capacity",
                                    levels=levels, cmap="Reds",
                                    vmin=vmin, vmax=vmax, norm=LogNorm())

    cbar.set_label("heat capacity")
    # `norm=LogNorm()` above sets the formatter to a LogFormatter by default
    cbar.formatter = ScalarFormatter()

    # Adjust size of plot to viewport to prevent clipping of the legend.
    plt.tight_layout()
    plt.show()


def read_true_longitude(filepath):
    return np.genfromtxt(filepath, dtype=np.float64)


def insolation(latitude, true_longitude, solar_constant, eccentricity,
               obliquity, precession_distance):
    # Determine if there is no sunset or no sunrise.
    sin_delta = np.sin(obliquity) * np.sin(true_longitude)
    cos_delta = np.sqrt(1 - sin_delta ** 2)
    tan_delta = sin_delta / cos_delta

    # Note that z can be +-infinity.
    # This is not a problem, as it is only used for the comparison with +-1.
    # We will never enter the `else` case below if z is +-infinity.
    z = -np.tan(latitude) * tan_delta

    if z >= 1:
        # Latitude where there is no sunrise
        return 0.0
    else:
        rho = ((1 - eccentricity * np.cos(true_longitude - precession_distance))
               / (1 - eccentricity ** 2)) ** 2

        if z <= -1:
            # Latitude where there is no sunset
            return solar_constant * rho * np.sin(latitude) * sin_delta
        else:
            h0 = np.arccos(z)
            second_term = h0 * np.sin(latitude) * sin_delta + np.cos(latitude) * cos_delta * np.sin(h0)
            return solar_constant * rho / np.pi * second_term


def calc_solar_forcing(albedo, true_longitudes, solar_constant=1371.685,
                       eccentricity=0.01674, obliquity=0.409253,
                       precession_distance=1.783037):
    def solar_forcing(theta, true_longitude, albedo_loc):
        s = insolation(theta, true_longitude, solar_constant, eccentricity,
                       obliquity, precession_distance)
        a_c = 1 - albedo_loc

        return s * a_c

    # Latitude values at the grid points
    nlatitude, nlongitude = albedo.shape
    y_lat = np.linspace(np.pi / 2, -np.pi / 2, nlatitude)

    return np.array([[[solar_forcing(y_lat[i], true_longitude, albedo[i, j])
                       for true_longitude in true_longitudes]
                      for j in range(nlongitude)]
                     for i in range(nlatitude)])


def plot_solar_forcing(solar_forcing, timestep, show_plot=False):
    vmin = np.amin(solar_forcing)
    vmax = np.amax(solar_forcing) * 1.05
    levels = np.linspace(vmin, vmax, 200)

    # Reuse plotting function from milestone 1.
    ntimesteps = solar_forcing.shape[2]
    day = (np.int_(timestep / ntimesteps * 365) + 80) % 365
    cbar = plot_robinson_projection(solar_forcing[:, :, timestep],
                                    f"Solar Forcing for Day {day}",
                                    levels=levels, cmap="gist_heat",
                                    vmin=vmin, vmax=vmax)
    cbar.set_label("solar forcing")

    # Adjust size of plot to viewport to prevent clipping of the legend.
    plt.tight_layout()

    filename = 'solar_forcing_{}.png'.format(timestep)
    plt.savefig(filename, dpi=300)
    if show_plot:
        plt.show()
    plt.close()

    return filename


# Run code
if __name__ == '__main__':
    geo_dat_ = read_geography("input/The_World128x65.dat")

    # Plot albedo
    albedo_ = calc_albedo(geo_dat_)
    plot_albedo(albedo_)

    # Plot heat capacity
    heat_capacity_ = calc_heat_capacity(geo_dat_)
    plot_heat_capacity(heat_capacity_)

    # Compute solar forcing
    true_longitude_ = read_true_longitude("input/True_Longitude.dat")
    solar_forcing_ = calc_solar_forcing(albedo_, true_longitude_)

    # Plot solar forcing for each time step
    filenames = []
    for ts in range(len(true_longitude_)):
        filename_ = plot_solar_forcing(solar_forcing_, ts, show_plot=False)
        filenames.append(filename_)

    # Build GIF
    frames = [imageio.v3.imread(filename_) for filename_ in filenames]
    imageio.mimsave("solar_forcing.gif", frames)

    # Remove files
    for filename_ in set(filenames):
        os.remove(filename_)
