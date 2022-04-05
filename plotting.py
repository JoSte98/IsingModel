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
    sigmas = []
    with open(name_of_file, "r") as file:
        lines = file.read().split('\n')
        del lines[-1]

        for line in lines:
            T, m, sigma = line.split(" ")
            temperatures.append(float(T))
            magnetizations.append(float(m))
            sigmas.append(float(sigma))

    fig, ax = plt.subplots()

    ax.errorbar(temperatures, magnetizations, yerr=sigmas, fmt='ro', label="Measurement")
    t_range_l = np.linspace(0.1, 2.269, 1000)
    t_range_h = np.linspace(2.269, temperatures[-1], 1000)
    ax.plot(t_range_l, (1 - t_range_l/2.269)**(1/8), lw=2.0, label=r'Phenomenological scaling with $\beta=1/8$',c="b")
    ax.plot(t_range_l, (1 - 1/(np.sinh(2/t_range_l)**4))**(1/8), lw=2.0, label="Exact Onsager solution",c="orange")
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

def plot_temp_susceptibility(name_of_file):
    """

        :param name_of_file:
        :return:
        """
    temperatures = []
    chis = []
    sigmas = []
    with open(name_of_file, "r") as file:
        lines = file.read().split('\n')
        del lines[-1]

        for line in lines:
            T, chi, sigma = line.split(" ")
            temperatures.append(float(T))
            chis.append(float(chi))
            sigmas.append(float(sigma))

    fig, ax = plt.subplots()

    ax.errorbar(temperatures, sigmas, yerr=sigmas, fmt='ro', label="Measurement")

    ax.set_xlim((0.0, temperatures[-1] + 0.2))
    ax.set_xlabel(r"Temperature $\frac{k_B T}{J}$", fontsize=18)
    ax.set_ylabel(r"Magnetic susceptibility $\chi_M$", fontsize=18)

    ax.legend(loc="best", fontsize=16)
    ax.grid(visible=True)

    plt.show()

def plot_temp_specific_heat(name_of_file):
    """

        :param name_of_file:
        :return:
        """
    temperatures = []
    Cs = []
    sigmas = []
    with open(name_of_file, "r") as file:
        lines = file.read().split('\n')
        del lines[-1]

        for line in lines:
            T, C, sigma = line.split(" ")
            temperatures.append(float(T))
            Cs.append(float(C))
            sigmas.append(float(sigma))

    fig, ax = plt.subplots()

    ax.errorbar(temperatures, Cs, yerr=sigmas, fmt='ro', label="Measurement")

    ax.set_xlim((0.0, temperatures[-1] + 0.2))
    ax.set_xlabel(r"Temperature $\frac{k_B T}{J}$", fontsize=18)
    ax.set_ylabel(r"Specific heat $C$", fontsize=18)

    ax.legend(loc="best", fontsize=16)
    ax.grid(visible=True)

    plt.show()

def plot_temp_energy(name_of_file):
    """

        :param name_of_file:
        :return:
        """
    temperatures = []
    es = []
    sigmas = []
    with open(name_of_file, "r") as file:
        lines = file.read().split('\n')
        del lines[-1]

        for line in lines:
            T, e, sigma = line.split(" ")
            temperatures.append(float(T))
            es.append(float(e))
            sigmas.append(float(sigma))

    fig, ax = plt.subplots()

    ax.errorbar(temperatures, es, yerr=sigmas, fmt='ro', label="Measurement")

    ax.set_xlim((0.0, temperatures[-1] + 0.2))
    ax.set_xlabel(r"Temperature $\frac{k_B T}{J}$", fontsize=18)
    ax.set_ylabel(r"Energy per spin $e$", fontsize=18)

    ax.legend(loc="best", fontsize=16)
    ax.grid(visible=True)

    plt.show()

def save_temp_magn(T, m, sigma, name_of_file):
    """

    :param T:
    :param m:
    :return:
    """
    with open(name_of_file, "a") as file:
        file.write("%f %f %f\n" % (T, m, sigma))

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

def save_temp_susceptibility(T, chi, sigma, name_of_file):
    """

    :param T:
    :param chi:
    :return:
    """
    with open(name_of_file, "a") as file:
        file.write("%f %f %f\n" % (T, chi, sigma))

    return 0

def save_temp_specific_heat(T, C, sigma, name_of_file):
    """

    :param T:
    :param C:
    :return:
    """
    with open(name_of_file, "a") as file:
        file.write("%f %f %f\n" % (T, C, sigma))

    return 0

def save_temp_energy(T, e, sigma, name_of_file):
    """

    :param T:
    :param e:
    :return:
    """
    with open(name_of_file, "a") as file:
        file.write("%f %f %f\n" % (T, e, sigma))

    return 0