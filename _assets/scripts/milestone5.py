import matplotlib.pyplot as plt
import numpy as np
import scipy
from numba import njit

from milestone1 import read_geography, plot_robinson_projection
from milestone2 import calc_albedo, calc_heat_capacity, calc_solar_forcing, read_true_longitude
from milestone3 import calc_area, calc_mean, calc_radiative_cooling_co2


class Mesh:
    def __init__(self, geo_dat):
        self.n_latitude, self.n_longitude = geo_dat.shape
        self.ndof = self.n_latitude * self.n_longitude
        self.h = np.pi / (self.n_latitude - 1)

        self.area = calc_area(geo_dat)
        self.geom = np.sin(0.5 * self.h) / self.area[0]

        self.csc2 = np.array([1 / np.sin(self.h * j) ** 2 for j in range(1, self.n_latitude - 1)])
        self.cot = np.array([1 / np.tan(self.h * j) for j in range(1, self.n_latitude - 1)])


def calc_diffusion_coefficients(geo_dat):
    nlatitude, nlongitude = geo_dat.shape

    coeff_ocean_poles = 0.40
    coeff_ocean_equator = 0.65
    coeff_equator = 0.65
    coeff_north_pole = 0.28
    coeff_south_pole = 0.20

    def diffusion_coefficient(j, i):
        # Compute the j value of the equator
        j_equator = int(nlatitude / 2)

        theta = np.pi * j / np.real(nlatitude - 1)
        colat = np.sin(theta) ** 5

        geo = geo_dat[j, i]
        if geo == 5:  # ocean
            return coeff_ocean_poles + (coeff_ocean_equator - coeff_ocean_poles) * colat
        else:  # land, sea ice, etc
            if j <= j_equator:  # northern hemisphere
                # on the equator colat=1 -> coefficients for norhern/southern hemisphere cancels out
                return coeff_north_pole + (coeff_equator - coeff_north_pole) * colat
            else:  # southern hemisphere
                return coeff_south_pole + (coeff_equator - coeff_south_pole) * colat

    return np.array(
        [[diffusion_coefficient(j, i) for i in range(nlongitude)] for j in range(nlatitude)]
    )


def calc_diffusion_operator(mesh, diffusion_coeff, temperature):
    h = mesh.h
    area = mesh.area
    n_latitude = mesh.n_latitude
    n_longitude = mesh.n_longitude
    csc2 = mesh.csc2
    cot = mesh.cot

    return calc_diffusion_operator_inner(h, area, n_latitude, n_longitude, csc2, cot, diffusion_coeff, temperature)


# This function is slow. It takes about 22ms on my machine, which is a lot for such a simple function.
# The problem is that Python is an interpreted (as opposed to compiled) language and therefore very slow.
# We could try to outsource most of that code to Numpy functions to speed things up, but that is not easy.
# Instead, we will use Numba, which is a JIT (just-in-time) compiler for Python to compile this function to fast
# machine code.
# With Numba, we can just add `@njit` to the function to tell Numba to compile it directly to machine code.
# Note, however, that this does not work with classes like our `Mesh`.
# Therefore, we cannot pass `mesh` to this function. Instead, we have to unpack all fields of `Mesh` above and
# add a wrapper function that does the unpacking.
#
# Times on my machine:
# Without `@njit`: 22ms (185s total for 8320 calls when calculating the Jacobian).
# With `@njit`: 25µs (209ms total for 8320 calls when calculating the Jacobian).
#
# This gives us a ~1000x speedup! Just by slightly rewriting the function and adding a decorator.
# Of course, we could try to further optimize this function, but that would be a waste of time, as it's not
# performance-relevant anymore.
@njit
def calc_diffusion_operator_inner(h, area, n_latitude, n_longitude, csc2, cot, diffusion_coeff, temperature):
    result = np.zeros(diffusion_coeff.shape)

    # North Pole
    factor = np.sin(h / 2) / (4 * np.pi * area[0])
    result[0, :] = factor * 0.5 * np.dot(diffusion_coeff[0, :] + diffusion_coeff[1, :],
                                         temperature[1, :] - temperature[0, :])

    # South Pole
    factor = np.sin(h / 2) / (4 * np.pi * area[-1])
    result[-1, :] = factor * 0.5 * np.dot(diffusion_coeff[-1, :] + diffusion_coeff[-2, :],
                                          temperature[-2, :] - temperature[-1, :])

    for i in range(n_longitude):
        # Only loop over inner nodes
        for j in range(1, n_latitude - 1):
            # There are the special cases of i=0 and i=n_longitude-1.
            # We have a periodic boundary condition, so for i=0, we want i-1 to be the last entry.
            # This happens automatically in Python when i=-1.
            # For i=n_longitude-1, we want i+1 to be 0.
            # For this, we define a variable ip1 (i plus 1) to avoid duplicating code.
            if i == n_longitude - 1:
                ip1 = 0
            else:
                ip1 = i + 1

            # Note that csc2 does not contain the values at the poles
            factor1 = csc2[j - 1] / h ** 2
            term1 = factor1 * (-2 * diffusion_coeff[j, i] * temperature[j, i] +
                               (diffusion_coeff[j, i] - 0.25 * (diffusion_coeff[j, ip1] - diffusion_coeff[j, i - 1])) *
                               temperature[j, i - 1] +
                               (diffusion_coeff[j, i] + 0.25 * (diffusion_coeff[j, ip1] - diffusion_coeff[j, i - 1])) *
                               temperature[j, ip1])

            factor2 = 1 / h ** 2
            term2 = factor2 * (-2 * diffusion_coeff[j, i] * temperature[j, i] +
                               (diffusion_coeff[j, i] - 0.25 *
                                (diffusion_coeff[j + 1, i] - diffusion_coeff[j - 1, i])) *
                               temperature[j - 1, i]
                               + (diffusion_coeff[j, i] + 0.25 *
                                  (diffusion_coeff[j + 1, i] - diffusion_coeff[j - 1, i])) *
                               temperature[j + 1, i])

            term3 = cot[j - 1] * diffusion_coeff[j, i] * 0.5 / h * (temperature[j + 1, i] - temperature[j - 1, i])

            result[j, i] = term1 + term2 + term3

    return result


def calc_operator_ebm_2d(temperature, mesh, diffusion_coeff, heat_capacity):
    diffusion_op = calc_diffusion_operator(mesh, diffusion_coeff, temperature)

    return (diffusion_op - 2.15 * temperature) / heat_capacity


def calc_source_terms_ebm_2d(heat_capacity, solar_forcing, radiative_cooling):
    return (solar_forcing - radiative_cooling) / heat_capacity


def calc_rhs_ebm_2d(temperature, mesh, diffusion_coeff, heat_capacity, solar_forcing, radiative_cooling):
    operator = calc_operator_ebm_2d(temperature, mesh, diffusion_coeff, heat_capacity)
    source_terms = calc_source_terms_ebm_2d(heat_capacity, solar_forcing, radiative_cooling)

    return operator + source_terms


def timestep_euler_forward_2d(temperature, t, delta_t,
                              mesh, diffusion_coeff, heat_capacity, solar_forcing, radiative_cooling):
    # Note that this function modifies the first argument instead of returning the result
    temperature[:, :, t] = temperature[:, :, t - 1] + delta_t * calc_rhs_ebm_2d(temperature[:, :, t - 1], mesh,
                                                                                diffusion_coeff, heat_capacity,
                                                                                solar_forcing[:, :, t - 1],
                                                                                radiative_cooling)


def compute_equilibrium_2d(timestep_function, mesh, diffusion_coeff, heat_capacity, solar_forcing, radiative_cooling,
                           max_iterations=100, rel_error=2e-5, verbose=True):
    # Number of time steps per year
    ntimesteps = solar_forcing.shape[2]

    # Step size
    delta_t = 1 / ntimesteps

    # We start with a constant temperature of 0 in every grid point throughout the year
    temperature = np.zeros((mesh.n_latitude, mesh.n_longitude, ntimesteps))

    # Area-mean in every time step
    area_mean_temp = np.zeros(ntimesteps)

    # Average temperature over all time steps from the previous iteration to approximate the error
    old_avg_temperature = 0

    for i in range(max_iterations):
        for t in range(ntimesteps):
            timestep_function(temperature, t, delta_t,
                              mesh, diffusion_coeff, heat_capacity, solar_forcing, radiative_cooling)
            area_mean_temp[t] = calc_mean(temperature[:, :, t], mesh.area)

        avg_temperature = np.sum(area_mean_temp) / ntimesteps
        verbose and print(f"Average annual temperature in iteration {i} is {avg_temperature}")

        if np.abs(avg_temperature - old_avg_temperature) < rel_error:
            # We can assume that the error is sufficiently small now.
            verbose and print("Equilibrium reached!")
            break
        else:
            old_avg_temperature = avg_temperature

    return temperature, area_mean_temp


def plot_diffusion_coefficient(diffusion_coeff):
    vmin = np.amin(diffusion_coeff)
    vmax = np.amax(diffusion_coeff)

    levels = np.linspace(vmin, vmax, 100)

    # Reuse plotting function from milestone 1.
    cbar = plot_robinson_projection(diffusion_coeff, "Diffusion Coefficients of the 2D EBM",
                                    levels=levels, cmap="cividis",
                                    vmin=vmin, vmax=vmax)
    cbar.set_ticks([vmin, 0.5 * (vmin + vmax), vmax])
    cbar.set_label("diffusion coefficient")

    # Adjust size of plot to viewport to prevent clipping of the legend.
    plt.tight_layout()
    plt.show()


def calc_jacobian_ebm_2d(mesh, diffusion_coeff, heat_capacity):
    jacobian = np.zeros((mesh.ndof, mesh.ndof))
    test_temperature = np.zeros(diffusion_coeff.shape)

    index = 0
    for j in range(mesh.n_latitude):
        for i in range(mesh.n_longitude):
            test_temperature[j, i] = 1.0
            op = calc_operator_ebm_2d(test_temperature, mesh, diffusion_coeff, heat_capacity)

            # Convert matrix to vector.
            # Note that this must be compatible with the loop order, so that `index` is correct.
            # `flatten` works row-wise, so we must loop over rows first in order to get the correct Jacobian.
            # To be precise, `test_temperature.flatten()` must be the `index`-th unit vector.
            jacobian[:, index] = op.flatten()

            # Reset test_temperature
            test_temperature[j, i] = 0.0
            index += 1

    return jacobian


# Run code
if __name__ == '__main__':
    geo_dat_ = read_geography("input/The_World128x65.dat")
    mesh_ = Mesh(geo_dat_)

    albedo = calc_albedo(geo_dat_)
    heat_capacity_ = calc_heat_capacity(geo_dat_)

    # Compute solar forcing
    true_longitude = read_true_longitude("input/True_Longitude.dat")
    solar_forcing_ = calc_solar_forcing(albedo, true_longitude)

    # Compute and plot diffusion coefficient
    diffusion_coeff_ = calc_diffusion_coefficients(geo_dat_)
    plot_diffusion_coefficient(diffusion_coeff_)

    co2_ppm = 315.0
    radiative_cooling_ = calc_radiative_cooling_co2(co2_ppm)

    compute_equilibrium_2d(timestep_euler_forward_2d, mesh_, diffusion_coeff_, heat_capacity_, solar_forcing_,
                           radiative_cooling_, max_iterations=2)

    # The Jacobian has three diagonals of non-zero entries and two blocks of non-zero entries for the poles.
    # We only show a subset of the entries, otherwise the two side diagonals are not visible.
    jacobian_ = calc_jacobian_ebm_2d(mesh_, diffusion_coeff_, heat_capacity_)
    plt.spy(jacobian_[0:300, 0:300])
    plt.show()

    # Use Scipy to efficiently compute the eigenvalues of this sparse matrix.
    # Converting the matrix to a sparse format first makes `eigs` ~100x faster.
    eigenvalues = scipy.sparse.linalg.eigs(scipy.sparse.csr_matrix(jacobian_), return_eigenvectors=False)

    # The maximum absolute value of the eigenvalues is too high for an efficient explicit time integration scheme.
    print(np.real(max(max(eigenvalues), -min(eigenvalues))))
