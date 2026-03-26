outdir = joinpath(@__DIR__, "..", "website", "_assets", "milestone_results")
mkpath(outdir)

include("milestone1.jl")
plot_geo_ = milestone1()
savefig(plot_geo_, joinpath(outdir, "milestone1_geo.png"))

include("milestone2.jl")
# Test `milestone2()` function and ignore the output
plot_albedo_, plot_heat_capacity_, gif_solar_forcing = milestone2()

savefig(plot_albedo_, joinpath(outdir, "milestone2_albedo.png"))
savefig(plot_heat_capacity_, joinpath(outdir, "milestone2_heat_capacity.png"))

# Move GIF file to output directory
mv(gif_solar_forcing.filename, joinpath(outdir, "milestone2_solar_forcing.gif"), force=true)

include("milestone3.jl")
plot_forward, plot_backward = milestone3()

savefig(plot_forward, joinpath(outdir, "milestone3_forward.png"))

include("milestone4.jl")
plot_mean, plot_pointwise, plot_cologne, gif_annual_temperature = milestone4()

savefig(plot_mean, joinpath(outdir, "milestone4_mean.png"))
savefig(plot_pointwise, joinpath(outdir, "milestone4_pointwise.png"))
savefig(plot_cologne, joinpath(outdir, "milestone4_cologne.png"))

# Move GIF file to output directory
mv(gif_annual_temperature.filename, joinpath(outdir, "milestone4_annual_temperature.gif"),
   force=true)

include("milestone5.jl")
plot_diffusion_coeff, plot_jacobian = milestone5()
savefig(plot_diffusion_coeff, joinpath(outdir, "milestone5_diffusion_coeff.png"))
savefig(plot_jacobian, joinpath(outdir, "milestone5_jacobian.png"))

include("milestone6.jl")
(plot_mean_temperature, plot_temperature_, plot_cologne, plot_temperature_co2,
 gif_temperature) = milestone6()
savefig(plot_temperature_, joinpath(outdir, "milestone6_temperature.png"))
savefig(plot_cologne, joinpath(outdir, "milestone6_cologne.png"))
savefig(plot_temperature_co2, joinpath(outdir, "milestone6_temperature_co2.png"))

# Move GIF file to output directory
mv(gif_temperature.filename, joinpath(outdir, "milestone6_annual_temperature.gif"),
   force=true)
