# Solve the 2D EBM

# Import everything
using LinearAlgebra: lu, ldiv!, UniformScaling
using SparseArrays: sparse
include(joinpath(@__DIR__,"..","milestone-1","geo.jl"))
include(joinpath(@__DIR__,"..","milestone-2","solarforcing.jl"))
include(joinpath(@__DIR__,"..","milestone-3","onedimensionalEBM.jl"))
include(joinpath(@__DIR__,"..","milestone-4","mesh.jl"))
include(joinpath(@__DIR__,"..","milestone-4","diffusioncoefficients.jl"))
include(joinpath(@__DIR__,"..","milestone-5","numerical_jacobian.jl"))


function milestone6()
    # Define parameters
    num_steps_year = 48
    dt = 1.0/num_steps_year
    tol = 2.0e-5

    ## Bla
    radiative_cooling_co2 = calc_radiative_cooling_co2(315.0)

    mesh = Mesh()
    radiative_cooling_feedback = 2.15   #[W/m^2/°C]: sensitivity of the seasonal cycle and annual change in the forcing agents

    geography = read_geography()
    diffusion_coeff = calc_diffusion_coefficients(geography,mesh.nx,mesh.ny)
    heatcap = calc_heat_capacity(geography)

    albedo = calc_albedo(geography,mesh.nx,mesh.ny)
    
    co_albedo = 1.0 .- albedo

    solar_forcing = calc_solar_forcing(co_albedo,mesh.nx,mesh.ny, num_steps_year)


    area = calc_area()

    # Solve the 2D EBM using forward Euler
    #####################################
    println("Computing with backward Euler")

    average_temperature = 0.0
    temperature = zeros(Float64,mesh.nx,mesh.ny,num_steps_year)  
    temperature.=5.0
    operator = zeros(Float64,mesh.nx,mesh.ny)  
    for year=1:100
        average_temperature_old = average_temperature
        average_temperature = 0.0
        
        for time_step=1:num_steps_year

            old_time_step = (time_step == 1) ? num_steps_year : time_step - 1

            compute_EBM_operator!(operator,mesh,radiative_cooling_feedback,diffusion_coeff,temperature[:,:,old_time_step],heatcap)    
            operator+= compute_EBM_source(solar_forcing[:,:,old_time_step],radiative_cooling_co2,heatcap)

            temperature[:,:,time_step] = temperature[:,:,old_time_step] + dt * operator

            monthly_average_temp = calc_mean(temperature[:,:,time_step], area)

            println(time_step," ",monthly_average_temp)

            average_temperature += monthly_average_temp
        end

        average_temperature *= dt

        println(average_temperature)
        if abs(average_temperature-average_temperature_old) < tol || isnan(average_temperature)
            break
        end
    end

    println("Computing with backward Euler")
    # Solve the 2D EBM using backward Euler
    #######################################

    # Compute jacobian
    jacobian = computeJacobian(mesh,radiative_cooling_feedback,diffusion_coeff,heatcap)

    jac_sparse = sparse(jacobian)
    lu_decomposition = lu(UniformScaling(1.0)-dt*jac_sparse)

    # Initializations
    average_temperature = 0.0
    temperature = zeros(Float64,mesh.nx,mesh.ny,num_steps_year)  
    temperature.=5.0

    temperature_vec = zeros(Float64,mesh.dof) 

    for year=1:100
        average_temperature_old = average_temperature
        average_temperature = 0.0
        
        for time_step=1:num_steps_year

            old_time_step = (time_step == 1) ? num_steps_year : time_step - 1

            rhs = temperature[:,:,old_time_step]
            rhs+= dt*compute_EBM_source(solar_forcing[:,:,old_time_step],radiative_cooling_co2,heatcap)

            #temperature[:,:,time_step] = temperature[:,:,old_time_step] + dt * operator

            ldiv!(temperature_vec, lu_decomposition, vec(rhs))
            
            temperature[:,:,time_step] = reshape(temperature_vec,(128,65))

            monthly_average_temp = calc_mean(temperature[:,:,time_step], area)

            #println(time_step," ",monthly_average_temp)

            average_temperature += monthly_average_temp
        end

        average_temperature *= dt
        println(year, " ", average_temperature)

        if abs(average_temperature-average_temperature_old) < tol
            break
        end
    end

    # Solve the 2D EBM using Crank-Nicolson
    #######################################
    println("Computing with Crank-Nicolson")

    # Compute jacobian
    jacobian = computeJacobian(mesh,radiative_cooling_feedback,diffusion_coeff,heatcap)

    jac_sparse = sparse(jacobian)
    lu_decomposition = lu(UniformScaling(1.0)-0.5*dt*jac_sparse)

    # Initializations
    average_temperature = 0.0
    temperature = zeros(Float64,mesh.nx,mesh.ny,num_steps_year)  
    temperature.=5.0

    temperature_vec = zeros(Float64,mesh.dof) 
    rhs = zeros(Float64,mesh.nx,mesh.ny)

    println(0, " ", calc_mean(temperature[:,:,num_steps_year], area))
    for year=1:100
        average_temperature_old = average_temperature
        average_temperature = 0.0
        
        for time_step=1:num_steps_year

            old_time_step = (time_step == 1) ? num_steps_year : time_step - 1
            
            compute_EBM_operator!(rhs,mesh,radiative_cooling_feedback,diffusion_coeff,temperature[:,:,old_time_step],heatcap)
            rhs=0.5*dt*rhs + temperature[:,:,old_time_step]
            rhs+= 0.5*dt*(compute_EBM_source(solar_forcing[:,:,old_time_step],radiative_cooling_co2,heatcap) + compute_EBM_source(solar_forcing[:,:,time_step],radiative_cooling_co2,heatcap))

            #temperature[:,:,time_step] = temperature[:,:,old_time_step] + dt * operator

            ldiv!(temperature_vec, lu_decomposition, vec(rhs))
            
            temperature[:,:,time_step] = reshape(temperature_vec,(128,65))

            monthly_average_temp = calc_mean(temperature[:,:,time_step], area)

            #println(time_step," ",monthly_average_temp)

            average_temperature += monthly_average_temp
        end

        average_temperature *= dt
        println(year, " ", average_temperature)

        if abs(average_temperature-average_temperature_old) < tol
            break
        end
    end

    return temperature, average_temperature, solar_forcing, albedo
end

function computeJacobian(mesh,radiative_cooling_feedback,diffusion_coeff,heatcap)
    jacobian = zeros(Float64,mesh.nx*mesh.ny,mesh.nx*mesh.ny)
    Fk = zeros(Float64,mesh.nx,mesh.ny)
    temperature_jac = zeros(Float64,mesh.nx,mesh.ny)  
    k=0
    for j=1:mesh.ny
        for i=1:mesh.nx
            k+=1
            temperature_jac[i,j]=1.0
            compute_EBM_operator!(Fk,mesh,radiative_cooling_feedback,diffusion_coeff,temperature_jac,heatcap)    
            jacobian[:,k] = vec(Fk)
            temperature_jac[i,j]=0.0
        end
    end

    return jacobian
end

temperature, average_temperature, solar_forcing, albedo= milestone6()