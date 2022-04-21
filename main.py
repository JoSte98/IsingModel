"""
Main file.
@author: Magdalena and Johannes
"""
import numpy as np
from plotting import *
from IsingModel import IsingModel

##### Names of files where statictical values should be stored after simulation.

NAME_PHASE_FILE = "phase_points.txt"
NAME_CORR_FILE = "correlation_time.txt"
NAME_ENERGY_FILE = "energy.txt"
NAME_SPECIFIC_HEAT_FILE = "specific_heat.txt"
NAME_SUSCEPTIBILITY_FILE = "susceptibility.txt"

def main():
    """

    :return: 0 if successful
    """

##### Calculating the statistic values of the system of 50x50 ('lenght' in line 30) spins for temperatures T in range (1.0,4.0). External magnetic field is preliminarily set to 0 here. In order to change that, please change 'external_field' in line 30).
##### The minimal number of measurement blocks is set to 10 here. Tn order to change that, look at second parameter in the line 32.
##### Here by changing the 'external_field' you can also verify, that all statistical values are independent of the sign of the external magnetic field. In case you would like to see how does the system prefers magnetization aligned with the external magnetic field, uncomment the section bellow (lines 50-64) and try for different values of the external field.
    
    for T in np.arange(1, 4.01, 0.2):
        T = np.round(T, 2)
        model = IsingModel(temperature=T, length=50, init_type='close_to_opt', external_field=0.0)

        model.measure_stat_prop(800, 10, plot=False, save_magnetizations=False)
        save_temp_corr_time(T, model.tau, NAME_CORR_FILE)
        save_temp_property(T, model.final_magnetization, model.final_magnetization_sigma, NAME_PHASE_FILE)
        save_temp_property(T, model.final_energy, model.final_energy_sigma, NAME_ENERGY_FILE)
        save_temp_property(T, model.final_susceptibility, model.final_susceptibility_sigma, NAME_SUSCEPTIBILITY_FILE)
        save_temp_property(T, model.final_specific_heat, model.final_specific_heat_sigma, NAME_SPECIFIC_HEAT_FILE)
        
##### Plotting of all collected data for given temperature range

    plot_temp_corr_time(NAME_CORR_FILE)
    plot_phase_diagram(NAME_PHASE_FILE)
    plot_temp_energy(NAME_ENERGY_FILE)
    plot_temp_specific_heat(NAME_SPECIFIC_HEAT_FILE)
    plot_temp_susceptibility(NAME_SUSCEPTIBILITY_FILE)
    
#####--------------------------------------------------------------------------------------------------  
##### Uncomment this whole section in case you want to try plotting of state configurations for different temperatures. By changing the 'external_field' you can study which magnetization is prefered for which sign of the external magnetic field.
#####
####    states=[]
####    temper = np.arange(1.0, 4.01, 0.2)
####    for T in temper:
####        T = np.round(T, 2)
####        model = IsingModel(temperature=T, length=50, init_type='close_to_opt', external_field=0.0)
####        model.simulate(100, plot=False, save=False)
####        states.append(model.state)
####    fig, ax = plt.subplots(4, 4, sharex='col', sharey='row')
####    plt.suptitle('Configurations of the system for different temperatures')
####    for i in range(4):
####        for j in range(4):
####            im = ax[i, j].imshow(states[4*i+j])
####            ax[i,j].set_title(r'$\tilde{T}$=' + str(np.round(temper[4*i+j],2)))
####    plt.tight_layout()
####    plt.show()
#####----------------------------------------------------------------------------------------------------------------  
##### This section shows an example how to just simulate the system independently for 500 steps and then automaticaly plot and save to file the evolution of the magnetization and energy of the system after the equilibration phase.
#####
####    model = IsingModel(temperature=1.5, length=50, init_type='up', external_field=0.0)
####    model.simulate(500,plot=True, save=True)
#####-----------------------------------------------------------------------------------------------------------------
    return 0

    



if __name__ == "__main__":
    main()