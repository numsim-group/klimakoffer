outdir = joinpath(@__DIR__, "..", "website", "_assets", "milestone_results")
mkpath(outdir)

include("milestone1.jl")
plot_geo_ = milestone1()
savefig(plot_geo_, joinpath(outdir, "milestone1_geo.png"))

include("milestone2.jl")
# Test `milestone2()` function and ignore the output
plot_albedo_, plot_heat_capacity_, gif_solar_forcing = milestone2()

# Plot solar forcing for day 80 of the year
geo_dat = read_geography(joinpath(@__DIR__, "input", "The_World128x65.dat"))

# Plot albedo
albedo = calc_albedo(geo_dat)

# Compute solar forcing
true_longitude = read_true_longitude(joinpath(@__DIR__, "input", "True_Longitude.dat"))
solar_forcing = calc_solar_forcing(albedo, true_longitude)

plot_solar_forcing_80 = plot_solar_forcing(solar_forcing, 1)

savefig(plot_albedo_, joinpath(outdir, "milestone2_albedo.png"))
savefig(plot_heat_capacity_, joinpath(outdir, "milestone2_heat_capacity.png"))
savefig(plot_solar_forcing_80, joinpath(outdir, "milestone2_solar_forcing_day_80.png"))

# Move GIF file to output directory
mv(gif_solar_forcing.filename, joinpath(outdir, "milestone2_solar_forcing.gif"), force=true)

include("milestone3.jl")
plot_forward, plot_backward = milestone3()

savefig(plot_forward, joinpath(outdir, "milestone3_forward.png"))

include("milestone4.jl")
(plot_mean, plot_pointwise, plot_cologne, plot_temperature_day_80,
 gif_annual_temperature) = milestone4()

savefig(plot_mean, joinpath(outdir, "milestone4_mean.png"))
savefig(plot_pointwise, joinpath(outdir, "milestone4_pointwise.png"))
savefig(plot_cologne, joinpath(outdir, "milestone4_cologne.png"))
savefig(plot_temperature_day_80, joinpath(outdir, "milestone4_temperature_day_80.png"))

# Move GIF file to output directory
mv(gif_annual_temperature.filename, joinpath(outdir, "milestone4_annual_temperature.gif"),
   force=true)

include("milestone5.jl")
plot_diffusion_coeff, plot_jacobian = milestone5()
savefig(plot_diffusion_coeff, joinpath(outdir, "milestone5_diffusion_coeff.png"))
savefig(plot_jacobian, joinpath(outdir, "milestone5_jacobian.png"))

# The Jacobian looks different in Python due to the row-major memory layout,
# so we also plot the Jacobian in row-major order to show how it looks in Python.
function calc_jacobian_row_major(mesh, diffusion_coeff, heat_capacity,
                                 radiative_cooling_feedback=2.15)
    jacobian = zeros((mesh.ndof, mesh.ndof))
    test_temperature = zeros(size(diffusion_coeff))
    op = similar(diffusion_coeff)

    index = 1
    # Here, the loop order is switched compared to `calc_jacobian_ebm_2d`
    for j in 1:(mesh.n_latitude), i in 1:(mesh.n_longitude)
        test_temperature[j, i] = 1.0
        # Use inplace version to avoid a lot of allocations
        calc_operator_ebm_2d!(op, test_temperature, mesh, diffusion_coeff, heat_capacity,
                              radiative_cooling_feedback)

        # Transpose before flattening to obtain row-major order
        jacobian[:, index] = vec(op')

        # Reset test_temperature
        test_temperature[j, i] = 0.0
        index += 1
    end

    return jacobian
end

function plot_jacobian_row_major()
    geo_dat = read_geography(joinpath(@__DIR__, "input", "The_World128x65.dat"))
    mesh = Mesh(geo_dat)

    albedo = calc_albedo(geo_dat)
    heat_capacity = calc_heat_capacity(geo_dat)

    # Compute solar forcing
    true_longitude = read_true_longitude(joinpath(@__DIR__, "input", "True_Longitude.dat"))
    solar_forcing = calc_solar_forcing(albedo, true_longitude)

    # Compute and plot diffusion coefficient
    diffusion_coeff = calc_diffusion_coefficients(geo_dat)

    co2_ppm = 315.0
    radiative_cooling = calc_radiative_cooling_co2(co2_ppm)

    compute_equilibrium_2d(timestep_euler_forward_2d, mesh, diffusion_coeff,
                           heat_capacity, solar_forcing,
                           radiative_cooling, max_iterations=2)

    # The Jacobian has three diagonals of non-zero entries and two blocks of non-zero entries for the poles.
    # We only show a subset of the entries, otherwise the two side diagonals are not visible.
    jacobian = calc_jacobian_row_major(mesh, diffusion_coeff, heat_capacity)

    # Plotting this matrix with a dense block is very slow with the PythonPlot backend,
    # so we switch to the GR backend for this plot.
    gr()
    return spy(jacobian[1:300, 1:300])
end
plot_jacobian_row_major_ = plot_jacobian_row_major()
savefig(plot_jacobian_row_major_, joinpath(outdir, "milestone5_jacobian_row_major.png"))

include("milestone6.jl")
(plot_mean_temperature, plot_temperature_, plot_cologne, plot_temperature_co2,
 plot_ziegler, plot_temperature_day_80, gif_temperature) = milestone6()
savefig(plot_temperature_, joinpath(outdir, "milestone6_temperature.png"))
savefig(plot_cologne, joinpath(outdir, "milestone6_cologne.png"))
savefig(plot_temperature_co2, joinpath(outdir, "milestone6_temperature_co2.png"))
savefig(plot_ziegler, joinpath(outdir, "milestone6_ziegler.png"))
savefig(plot_temperature_day_80, joinpath(outdir, "milestone6_temperature_day_80.png"))

# Move GIF file to output directory
mv(gif_temperature.filename, joinpath(outdir, "milestone6_annual_temperature.gif"),
   force=true)
