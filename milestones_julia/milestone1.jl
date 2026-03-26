using DelimitedFiles
using Plots
pythonplot()

# Creates an nlatitude x nlongitude = 65 x 128 array with integer digits decoding the geography.
function read_geography(filepath)
    return readdlm(filepath)
end

function robinson_projection(nlatitude, nlongitude)
    function x_fun(lon, lat)
        return lon / pi * (0.0379 * lat^6 - 0.15 * lat^4 - 0.367 * lat^2 + 2.666)
    end

    function y_fun(_, lat)
        return 0.96047 * lat - 0.00857 * sign(lat) * abs(lat)^6.41
    end

    # Longitude goes from -pi to pi (not included), latitude from pi/2 to -pi/2.
    # Remove endpoint of longitude to avoid overlap.
    x_lon = LinRange(-pi, pi, nlongitude + 1)[1:(end - 1)]

    # Latitude goes backwards because our plot starts at the bottom,
    # but the matrix is read from the top.
    y_lat = LinRange(pi / 2, -pi / 2, nlatitude)

    x = [x_fun(lon, lat) for lat in y_lat, lon in x_lon]
    y = [y_fun(lon, lat) for lat in y_lat, lon in x_lon]

    return x, y
end

function plot_geo(geo_dat)
    nlatitude, nlongitude = size(geo_dat)
    x, y = robinson_projection(nlatitude, nlongitude)

    # Unfortunately, this is a bit buggy. When passing a `cgrad` with `categorical=true`,
    # one would expect to only get the colors of this discrete `cgrad`, but for some reason,
    # that is not the case.
    # Instead, we got the colors that we want by experimenting with the levels.
    # We tried to make the range for ocean only slightly larger than the others to avoid
    # a weird looking colorbar.
    plot = contourf(x, y, geo_dat,
                    levels=[0.5, 1.7, 2.9, 4.1, 5.5],
                    clims=(1, 5),
                    aspect_ratio=1,
                    title="Earth Geography",
                    c=cgrad([:darkgreen, :lightsteelblue, :lavender, :navy]),
                    colorbar_ticks=([1.1, 2.3, 3.5, 4.8],
                                    ["land", "sea ice", "snow cover", "ocean"]),
                    axis=([], false),
                    dpi=300)

    return plot
end

# Run code
function milestone1()
    geo_dat = read_geography(joinpath(@__DIR__, "input", "The_World128x65.dat"))
    plot_geo(geo_dat)
end
