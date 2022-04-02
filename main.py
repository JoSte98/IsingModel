"""
@author: Magdalena and Johannes
"""
import numpy as np

from plotting import *
from IsingModel import IsingModel

NAME_PHASE_FILE = "phase_points.txt"
NAME_CORR_FILE = "correlation_time.txt"

def main():
    """

    :return: 0 if successful
    """

    for T in np.arange(1, 4.01, 0.2):
        T = np.round(T, 2)
        if T < 2.1:
            model = IsingModel(temperature=T, length=50, init_type='up')
        else:
            model = IsingModel(temperature=T, length=50, init_type='random')

        model.measure_corr_time(800, plot=False)
        save_temp_corr_time(T, model.tau, NAME_CORR_FILE)
        save_temp_magn(T, np.mean(model.magnetizations), NAME_PHASE_FILE)


    plot_temp_corr_time(NAME_CORR_FILE)
    plot_phase_diagram(NAME_PHASE_FILE)


    return 0



if __name__ == "__main__":
    main()