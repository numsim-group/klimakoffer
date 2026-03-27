import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from milestone3 import calc_rhs_system_backward_euler, calc_matrix_system_backward_euler, compute_equilibrium_EBM_0D

def compute_evolution_EBM_0D(mean_heat_capacity, mean_solar_forcing):
    # read co2 nasa data
    filepath_co2 = "input/co2_nasa.dat"
    nasa_dat = np.genfromtxt(filepath_co2)
    n_month,n_data = nasa_dat.shape
    years = int(n_month / 12) # assumption that only full years are available
    ntimesteps = len(mean_solar_forcing)
    n_average_temperatures = years
    n_annual_temperatures = years * ntimesteps
    first_year = int(nasa_dat[0,0])
    last_year = int(nasa_dat[n_month-1,0])
    
    # allocate
    average_temperatures = np.zeros(n_average_temperatures,dtype=np.float64)
    annual_temperatures = np.zeros(n_annual_temperatures,dtype=np.float64)
    co2_data = np.zeros(years,dtype=np.float64)
    
    #compute average of co2 per year
    for yr in range(years):
        for month in range(12):
            index = 12 * yr + month
            co2_data[yr] += nasa_dat[index,3]
        co2_data[yr] = co2_data[yr] / 12       
        start_index = 0 + ntimesteps * yr
        end_index = 47 + ntimesteps * yr
        average_temperature, annual_temperature = compute_equilibrium_EBM_0D(mean_heat_capacity, mean_solar_forcing, calc_matrix_system_forward_euler, calc_rhs_system_forward_euler, max_iterations=100, rel_error=1e-2, verbose=False, co2_ppm = co2_data[yr])
        average_temperatures[yr] = average_temperature
        annual_temperatures[start_index:end_index+1] = annual_temperature
    
    print('EVOLUTION FINISHED!')
    
    return average_temperatures, annual_temperatures, first_year, last_year

def compute_transient_EBM_0D(mean_heat_capacity, mean_solar_forcing):
    # read co2 nasa data
    filepath_co2 = "input/co2_nasa.dat"
    nasa_dat = np.genfromtxt(filepath_co2)
    n_month,n_data = nasa_dat.shape
    #print(n_month,n_data)
    years = int(n_month / 12) # assumption that only full years are available
    #print(years)
    ntimesteps = len(mean_solar_forcing)
    # delta t in the year
    delta_t = 1.0 / np.real(ntimesteps)
    
    n_average_temperatures = years
    n_annual_temperatures = years * ntimesteps
    first_year = int(nasa_dat[0,0])
    last_year = int(nasa_dat[n_month-1,0])
    
    # allocate
    average_temperatures = np.zeros(n_average_temperatures,dtype=np.float64)
    annual_temperatures = np.zeros(n_annual_temperatures,dtype=np.float64)
    co2_data = np.zeros(years,dtype=np.float64)
    
    
    # compute equilibrium for the first year   
    # compute average co2 for the first year
    for month in range(12):
        index = month
        co2_data[0] += nasa_dat[index,3]
        
    co2_data[0] = co2_data[0] / 12
    
    # get linear system matrix (only scalar quantity in the 0D EBM!)
    matrix_system = calc_matrix_system_backward_euler(delta_t,mean_heat_capacity, co2_ppm = co2_data[0])

    average_temperature, annual_temperature = compute_equilibrium_EBM_0D(mean_heat_capacity, mean_solar_forcing, max_iterations=100, rel_error=1e-2, verbose=False, co2_ppm = co2_data[0])
    average_temperatures[0] = average_temperature
    annual_temperatures[0:48] = annual_temperature
        
    #compute average of co2 per year
    for yr in range(1,years):
        # compute yearly co2 average
        for month in range(12):
            index = 12 * yr + month
            co2_data[yr] += nasa_dat[index,3]
          
        co2_data[yr] = co2_data[yr] / 12       
        
        average_temperature = 0.0
        for ts in range(ntimesteps):
            ts_old = ts - 1 #note that we abuse the python regulation that -1 gives the last entry! 
            # integrate solar forcing
            # get rhs of linear system
            rhs_system = calc_rhs_system_backward_euler(delta_t,mean_heat_capacity,  mean_solar_forcing[ts_old], mean_solar_forcing[ts], annual_temperature[ts_old], co2_ppm)
            annual_temperature[ts] = rhs_system / matrix_system
            average_temperature += annual_temperature[ts]
        
        average_temperature = average_temperature / np.real(ntimesteps)
        # store in array       
        start_index = 0 + ntimesteps * yr
        end_index = 47 + ntimesteps * yr
        average_temperatures[yr] = average_temperature
        annual_temperatures[start_index:end_index+1] = annual_temperature
    
    print('EVOLUTION FINISHED!')
    
    return average_temperatures, annual_temperatures, first_year, last_year
