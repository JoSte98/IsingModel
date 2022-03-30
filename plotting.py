"""

"""

import numpy as np
import matplotlib.pyplot as plt

def plot_phase_diagram(name_of_file):
    """

    :param name_of_file:
    :return:
    """
    temperatures = []
    magnetizations = []
    with open(name_of_file, "r") as file:
        lines = file.read().split('\n')
        del lines[-1]

        for line in lines:
            T, m = line.split(" ")
            temperatures.append(float(T))
            magnetizations.append(float(m))

    fig, ax = plt.subplots()

    print(temperatures)
    print(magnetizations)

    ax.scatter(temperatures, magnetizations, linewidths=2.0, label="Measurement")
    t_range = np.linspace(temperatures[0], temperatures[-1], 1000)
    ax.plot(t_range, np.sqrt(1 - t_range/2.269), lw=2.0, label="Theory")

    ax.set_xlim((0.0, temperatures[-1] + 0.2))

    ax.legend(loc="best", fontsize=16)

    plt.show()


def save_temp_magn(T, m, name_of_file):
    """

    :param T:
    :param m:
    :return:
    """
    with open(name_of_file, "a") as file:
        file.write("%f %f\n" % (T, m))

    return 0

