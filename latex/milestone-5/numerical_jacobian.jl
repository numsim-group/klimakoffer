using SparseArrays: sparse
using UnPack: @unpack

function compute_EBM_operator!(F,mesh,radiative_cooling_feedback,diffusion_coeff,temperature,heatcap)
    @unpack nx,ny,dof,h,geom,csc2,cot,area = mesh
    
    # Initializations
    sh2 = 1 / h^2
    sh   = 1 / h
    
    # Inner DOFs
    for j=2:ny-1 # Loop over latitude

        i = 1 # (periodic BC at lambda = -180°)
        # Add solar forcing, radiative cooling, and feedback effects
        ############################################################
        F[i,j] = - radiative_cooling_feedback * temperature[i,j]

        # Add head conduction terms
        ###########################
        # Second derivative term in longitude
        F[i,j] += csc2[j] * sh2 * (- 2 * diffusion_coeff[i,j] * temperature[i,j] 
                                    + (diffusion_coeff[i,j] - 0.25 * (diffusion_coeff[i+1,j] - diffusion_coeff[nx,j])) * temperature[nx,j] 
                                    + (diffusion_coeff[i,j] + 0.25 * (diffusion_coeff[i+1,j] - diffusion_coeff[nx,j])) * temperature[i+1,j] )
        # Second derivative term in latitude
        F[i,j] += sh2 * (- 2 * diffusion_coeff[i,j] * temperature[i,j] 
                            + (diffusion_coeff[i,j] - 0.25 * (diffusion_coeff[i,j+1] - diffusion_coeff[i,j-1])) * temperature[i,j-1] 
                            + (diffusion_coeff[i,j] + 0.25 * (diffusion_coeff[i,j+1] - diffusion_coeff[i,j-1])) * temperature[i,j+1] )
        # First derivative term in latitude
        F[i,j] += cot[j] * diffusion_coeff[i,j] * 0.5 * sh * (temperature[i,j+1] - temperature[i,j-1])

        for i=2:nx-1     # Loop over longitude
            
            # Add solar forcing, radiative cooling, and feedback effects
            ############################################################
            F[i,j] = - radiative_cooling_feedback * temperature[i,j]

            # Add head conduction terms
            ###########################
            # Second derivative term in longitude
            F[i,j] += csc2[j] * sh2 * (- 2 * diffusion_coeff[i,j] * temperature[i,j] 
                                       + (diffusion_coeff[i,j] - 0.25 * (diffusion_coeff[i+1,j] - diffusion_coeff[i-1,j])) * temperature[i-1,j] 
                                       + (diffusion_coeff[i,j] + 0.25 * (diffusion_coeff[i+1,j] - diffusion_coeff[i-1,j])) * temperature[i+1,j] )
            # Second derivative term in latitude
            F[i,j] += sh2 * (- 2 * diffusion_coeff[i,j] * temperature[i,j] 
                             + (diffusion_coeff[i,j] - 0.25 * (diffusion_coeff[i,j+1] - diffusion_coeff[i,j-1])) * temperature[i,j-1] 
                             + (diffusion_coeff[i,j] + 0.25 * (diffusion_coeff[i,j+1] - diffusion_coeff[i,j-1])) * temperature[i,j+1] )
            # First derivative term in latitude
            F[i,j] += cot[j] * diffusion_coeff[i,j] * 0.5 * sh * (temperature[i,j+1] - temperature[i,j-1])
            
        end

        i = nx # (periodic BC at lambda = 180°)
        # Add solar forcing, radiative cooling, and feedback effects
        ############################################################
        F[i,j] = - radiative_cooling_feedback * temperature[i,j]

        # Add head conduction terms
        ###########################
        # Second derivative term in longitude
        F[i,j] += csc2[j] * sh2 * (- 2 * diffusion_coeff[i,j] * temperature[i,j] 
                                    + (diffusion_coeff[i,j] - 0.25 * (diffusion_coeff[1,j] - diffusion_coeff[i-1,j])) * temperature[i-1,j] 
                                    + (diffusion_coeff[i,j] + 0.25 * (diffusion_coeff[1,j] - diffusion_coeff[i-1,j])) * temperature[1,j] )
        # Second derivative term in latitude
        F[i,j] += sh2 * (- 2 * diffusion_coeff[i,j] * temperature[i,j] 
                            + (diffusion_coeff[i,j] - 0.25 * (diffusion_coeff[i,j+1] - diffusion_coeff[i,j-1])) * temperature[i,j-1] 
                            + (diffusion_coeff[i,j] + 0.25 * (diffusion_coeff[i,j+1] - diffusion_coeff[i,j-1])) * temperature[i,j+1] )
        # First derivative term in latitude
        F[i,j] += cot[j] * diffusion_coeff[i,j] * 0.5 * sh * (temperature[i,j+1] - temperature[i,j-1])
    end

    # DOFs at the poles

    # North pole
    # **********
    # Add solar forcing, radiative cooling, and feedback effects
    ############################################################
    F[:,1] = - radiative_cooling_feedback .* temperature[:,1]
    # Add head conduction terms
    ###########################
    for k=1:nx
        mean_diffusion_coeff = (diffusion_coeff[k,1] * area[1] + diffusion_coeff[k,2] * area[2]) / (area[1] + area[2])
        geom_mean_diffusion_coeff = geom * mean_diffusion_coeff
        #= # Slow version
        F[:,1] .+= geom_mean_diffusion_coeff * temperature[k,2]
        F[:,1] -= geom_mean_diffusion_coeff .* temperature[:,1] =#
        # Slow version mod
        #= @views F[:,1] .+= geom_mean_diffusion_coeff * temperature[k,2]
        @views F[:,1] -= geom_mean_diffusion_coeff .* temperature[:,1] =#
        # Fast version 
        for i=1:nx
            F[i,1] += geom_mean_diffusion_coeff * temperature[k,2]
            F[i,1] -= geom_mean_diffusion_coeff * temperature[i,1]
        end 
    end

    # South pole
    # **********
    # Add solar forcing, radiative cooling, and feedback effects
    ############################################################
    F[:,ny] = - radiative_cooling_feedback .* temperature[:,ny]
    # Add head conduction terms
    ###########################
    for k=1:nx
        mean_diffusion_coeff = (diffusion_coeff[k,ny] * area[ny] + diffusion_coeff[k,ny-1] * area[ny-1]) / (area[ny] + area[ny-1])
        geom_mean_diffusion_coeff = geom * mean_diffusion_coeff
        #= # Slow version
        F[:,ny] -= geom_mean_diffusion_coeff .* temperature[:,ny] 
        F[:,ny] .+= geom_mean_diffusion_coeff * temperature[k,ny-1] =#
        # Slow version mod
        #= @views F[:,ny] -= geom_mean_diffusion_coeff .* temperature[:,ny] 
        @views F[:,ny] .+= geom_mean_diffusion_coeff * temperature[k,ny-1] =#
        # Fast version
        for i=1:nx
            F[i,ny] -= geom_mean_diffusion_coeff * temperature[i,ny] 
            F[i,ny] += geom_mean_diffusion_coeff * temperature[k,ny-1]
        end
    end
    
    F .= F./heatcap
end

function compute_EBM_source(solar_forcing,radiative_cooling_co2,heatcap)
    return (solar_forcing .- radiative_cooling_co2)./heatcap
end

function milestone5()
    include(joinpath(@__DIR__,"..","milestone-4","mesh.jl"))
    include(joinpath(@__DIR__,"..","milestone-2","solarforcing.jl"))
    include(joinpath(@__DIR__,"..","milestone-4","diffusioncoefficients.jl"))
    include(joinpath(@__DIR__,"..","milestone-1","geo.jl"))

    mesh = Mesh()
    radiative_cooling_feedback = 2.15   #[W/m^2/°C]: sensitivity of the seasonal cycle and annual change in the forcing agents
    
    geography = read_geography(joinpath(@__DIR__,"input","The_World128x65.dat"),mesh.nx,mesh.ny)
    heatcap = calc_heat_capacity(geography)

    diffusion_coeff = calc_diffusion_coefficients(geography,mesh.nx,mesh.ny)

    jacobian = zeros(Float64,mesh.nx*mesh.ny,mesh.nx*mesh.ny)

    temperature = zeros(Float64,mesh.nx,mesh.ny)  
    Fk = zeros(Float64,mesh.nx,mesh.ny)
    k=0
    for j=1:mesh.ny
        for i=1:mesh.nx
            k+=1
            temperature[i,j]=1
            compute_EBM_operator!(Fk,mesh,radiative_cooling_feedback,diffusion_coeff,temperature,heatcap)    
            jacobian[:,k] = vec(Fk)
            temperature[i,j]=0
        end
    end

    return jacobian

end


