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

    ax.scatter(temperatures, magnetizations, linewidths=1.5, c="r", label="Measurement")
    t_range_l = np.linspace(temperatures[0], 2.269, 1000)
    t_range_h = np.linspace(2.269, temperatures[-1], 1000)
    ax.plot(t_range_l, np.sqrt(1 - t_range_l/2.269), lw=2.0, label="Theory")
    ax.plot(t_range_h, np.zeros(len(t_range_h)), lw=2.0, c="b")

    ax.set_xlim((0.0, temperatures[-1] + 0.2))
    ax.set_xlabel(r"Temperature $\frac{k_B T}{J}$", fontsize=18)
    ax.set_ylabel(r"Magnetization per spin $\langle M\rangle$", fontsize=18)
    ax.grid(visible=True)

    ax.legend(loc="best", fontsize=16)

    plt.show()

def plot_temp_corr_time(name_of_file):
    """

        :param name_of_file:
        :return:
        """
    temperatures = []
    taus = []
    with open(name_of_file, "r") as file:
        lines = file.read().split('\n')
        del lines[-1]

        for line in lines:
            T, tau = line.split(" ")
            temperatures.append(float(T))
            taus.append(float(tau))

    fig, ax = plt.subplots()

    ax.scatter(temperatures, taus, linewidths=1.5, c="r", label="Measurement")

    ax.set_xlim((0.0, temperatures[-1] + 0.2))
    ax.set_xlabel(r"Temperature $\frac{k_B T}{J}$", fontsize=18)
    ax.set_ylabel(r"Correlation time $\tau$", fontsize=18)

    ax.legend(loc="best", fontsize=16)
    ax.grid(visible=True)

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

def save_temp_corr_time(T, tau, name_of_file):
    """

    :param T:
    :param tau:
    :return:
    """
    with open(name_of_file, "a") as file:
        file.write("%f %f\n" % (T, tau))

    return 0