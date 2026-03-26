using Interpolations;
using Plots
pyplot()


"""
read_geography()
Read input data to get the surface distribution.
"""
function read_geography(filepath=joinpath(@__DIR__,"input","The_World128x65.dat"),nlongitude=128,nlatitude=65)
    result = zeros(Int8,nlongitude,nlatitude)
    open(filepath) do fh
        for lat = 1:nlatitude
            if eof(fh) break end
               result[:,lat] = parse.(Int8,split(strip(readline(fh) ),r""))
            end
        end
return result
end


"""
apply_robinson_projection()
The Robinson projection is a map network design
based on reference points that have been empirically determined.
"""
function apply_robinson_projection(data)

    lat = [0 , 5 , 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90] * pi/180.0
    X   = [1.0000, 0.9986, 0.9954, 0.99  , 0.9822, 0.973 , 0.96  , 0.9427, 0.9216, 0.8962, 0.8679, 0.835 , 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322]
    Y   = [0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000]

    itpl_X = LinearInterpolation(lat,X);
    itpl_Y = LinearInterpolation(lat,Y);

    calc_x = (radius,longitude,latitude) -> 0.8487 * radius .* itpl_X.(abs.(latitude)) .* longitude;
    calc_y = (radius,longitude,latitude) -> 1.3523 * radius .* itpl_Y.(abs.(latitude)) .* 0.2536 .* sign.(latitude)

    nx = size(data,1)
    ny = size(data,2)

    longitude = LinRange(-pi,pi,nx) .* ones(ny)'
    latitude  = ones(nx) .* LinRange(-pi/2,pi/2,ny)'

    x = calc_x(1.0,longitude,latitude);
    y = calc_y(1.0,longitude,latitude);
    z = reverse(reverse(data),dims=1)

    x = 180.0*(x ./ maximum(x));
    y =  90.0*(y ./ maximum(y));

    return x,y,z

end


function plot_geo(geography)
 
    data = geography    
    x,y,z = apply_robinson_projection(data)

    title = string("Earth's surface");

    z = reverse(reverse(data),dims=1);
    plot = contourf(x,y,z,
    clims=(1,5),
    levels=LinRange(1,5,500),
    aspect_ratio=1,
    title=title,
    xlabel="longitude [°]",
    ylabel="latitude [°]",
    colorbar_title="surface types",
    c=:PuBu,
    size=(1024,640))

    savefig(plot,string("earth_surface.png"))

 end

function milestone1()
    nlongitude = 128
    nlatitude = Int(nlongitude/2+1) # 65
    geography = read_geography(joinpath(@__DIR__,"input","The_World128x65.dat"),nlongitude,nlatitude)

    plot_geo(geography)
end
