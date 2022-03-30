"""
@author: Magdalena and Johannes
"""
import numpy as np

from plotting import plot_phase_diagram, save_temp_magn

NAME_PHASE_FILE = "phase_points.txt"

def main():
    """

    :return:
    """
    from IsingModel import IsingModel

    model = IsingModel(temperature=3, length=50)
    model.simulate(20)
    model.plot_energies()

#    for T in np.arange(1, 4, 0.2):
#        model = IsingModel(temperature=T, length=20, init_type='up')
#        model.simulate(10, plot=False)
#
#        save_temp_magn(T, model.magnetizations[-1], NAME_PHASE_FILE)
#
#
#    plot_phase_diagram(NAME_PHASE_FILE)


    return 0



if __name__ == "__main__":
    main()