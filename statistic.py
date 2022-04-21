"""
Functions to get the Ising model properties for different temperatures.

@author: Magdalena and Johannes
"""
from IsingModel import IsingModel
import numpy as np
from plotting import *

def Ising_model_statistic(T_start, T_stop, length, NAME_PHASE_FILE, NAME_ENERGY_FILE, NAME_SUSCEPTIBILITY_FILE,NAME_SPECIFIC_HEAT_FILE, external_field=0.0):
    for T in np.arange(1, 4.01, 0.2):
        T = np.round(T, 2)
        model = IsingModel(temperature=T, length=length, init_type='close_to_opt', external_field=external_field)

        model.measure_stat_prop(800, 10, plot=False, save_magnetizations=False)
        save_temp_corr_time(T, model.tau, NAME_CORR_FILE)
        save_temp_property(T, model.final_magnetization, model.final_magnetization_sigma, NAME_PHASE_FILE)
        save_temp_property(T, model.final_energy, model.final_energy_sigma, NAME_ENERGY_FILE)
        save_temp_property(T, model.final_susceptibility, model.final_susceptibility_sigma, NAME_SUSCEPTIBILITY_FILE)
        save_temp_property(T, model.final_specific_heat, model.final_specific_heat_sigma, NAME_SPECIFIC_HEAT_FILE)
        
    return 0

def Ising_model_configurations(length, external_field=0.0):
    states=[]
    temper = np.arange(1.0, 4.01, 0.2)
    for T in temper:
        T = np.round(T, 2)
        model = IsingModel(temperature=T, length=length, init_type='close_to_opt', external_field=external_field)
        model.simulate(100, plot=False, save=False)
        states.append(model.state)
    fig, ax = plt.subplots(4, 4, sharex='col', sharey='row')
    plt.suptitle('Configurations of the system for different temperatures')
    for i in range(4):
        for j in range(4):
            im = ax[i, j].imshow(states[4*i+j])
            ax[i,j].set_title(r'$\tilde{T}$=' + str(np.round(temper[4*i+j],2)))
    plt.tight_layout()
    plt.show()
    return 0