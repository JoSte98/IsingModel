"""
@author: Magdalena and Johannes
"""
import numpy as np

from plotting import *
from IsingModel import IsingModel

NAME_PHASE_FILE = "phase_points.txt"
NAME_CORR_FILE = "correlation_time.txt"
NAME_ENERGY_FILE = "energy.txt"
NAME_SPECIFIC_HEAT_FILE = "specific_heat.txt"
NAME_SUSCEPTIBILITY_FILE = "susceptibility.txt"

def main():
    """

    :return: 0 if successful
    """

    for T in np.arange(1, 4.01, 0.2):
        T = np.round(T, 2)
        model = IsingModel(temperature=T, length=50, init_type='close_to_opt')

        model.measure_stat_prop(800, 10, plot=False)

        save_temp_corr_time(T, model.tau, NAME_CORR_FILE)
        save_temp_magn(T, model.final_magnetization, model.final_magnetization_sigma, NAME_PHASE_FILE)
        save_temp_energy(T, model.final_energy, model.final_energy_sigma, NAME_ENERGY_FILE)
        save_temp_susceptibility(T, model.final_susceptibility, model.final_susceptibility_sigma, NAME_SUSCEPTIBILITY_FILE)
        save_temp_specific_heat(T, model.final_specific_heat, model.final_specific_heat_sigma, NAME_SPECIFIC_HEAT_FILE)

    plot_temp_corr_time(NAME_CORR_FILE)
    plot_phase_diagram(NAME_PHASE_FILE)
    plot_temp_energy(NAME_ENERGY_FILE)
    plot_temp_specific_heat(NAME_SPECIFIC_HEAT_FILE)
    plot_temp_susceptibility(NAME_SUSCEPTIBILITY_FILE)


    return 0



if __name__ == "__main__":
    main()