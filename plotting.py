"""
External plotting and saving to a file functions for simulated statistical properties.

@author: Magdalena and Johannes
"""

import numpy as np
import matplotlib.pyplot as plt

def load_property_from_file(name_of_file):
    """
    Loads temperatures, means and errors of a given property (magnetization/energy/susceptibility/specific heat) from a file.
    
    :param name_of_file: Name of the file from which values are loaded.
    :return: Lists of temperatures, mean values, errors
    """
    temperatures = []
    means = []
    sigmas = []
    with open(name_of_file, "r") as file:
        lines = file.read().split('\n')
        del lines[-1]

        for line in lines:
            T, mean, sigma = line.split(" ")
            temperatures.append(float(T))
            means.append(float(mean))
            sigmas.append(float(sigma))
    
    return temperatures, means, sigmas

def plot_phase_diagram(name_of_file):
    """
    Plots phase diagram (magnetizations for different temperatures) of the system from values saved in a file.

    :param name_of_file: Name of the file from which values are loaded.
    :return: 0 if successful
    """
    temperatures, magnetizations, sigmas = load_property_from_file(name_of_file)

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
    
    return 0

def plot_temp_corr_time(name_of_file):
    """
    Plots correlation time of the system for different temperatures from values saved in a file.

    :param name_of_file: Name of the file from which values are loaded.
    :return: 0 if successful
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
    
    return 0

def plot_temp_susceptibility(name_of_file):
    """
    Plots susceptibility of the system for different temperatures from values saved in a file.

    :param name_of_file: Name of the file from which values are loaded.
    :return: 0 if successful
    """
    temperatures, susceptibilities, sigmas = load_property_from_file(name_of_file)

    fig, ax = plt.subplots()

    ax.errorbar(temperatures, susceptibilities, yerr=sigmas, fmt='ro', label="Measurement")

    ax.set_xlim((0.0, temperatures[-1] + 0.2))
    ax.set_xlabel(r"Temperature $\frac{k_B T}{J}$", fontsize=18)
    ax.set_ylabel(r"Magnetic susceptibility $\chi_M$", fontsize=18)

    ax.legend(loc="best", fontsize=16)
    ax.grid(visible=True)

    plt.show()
    
    return 0

def plot_temp_specific_heat(name_of_file):
    """
    Plots specific heat of the system for different temperatures from values saved in a file.

    :param name_of_file: Name of the file from which values are loaded.
    :return: 0 if successful
    """
    temperatures, heats, sigmas = load_property_from_file(name_of_file)

    fig, ax = plt.subplots()

    ax.errorbar(temperatures, heats, yerr=sigmas, fmt='ro', label="Measurement")

    ax.set_xlim((0.0, temperatures[-1] + 0.2))
    ax.set_xlabel(r"Temperature $\frac{k_B T}{J}$", fontsize=18)
    ax.set_ylabel(r"Specific heat $C$", fontsize=18)

    ax.legend(loc="best", fontsize=16)
    ax.grid(visible=True)

    plt.show()
    
    return 0

def plot_temp_energy(name_of_file):
    """
    Plots energy per spin of the system for different temperatures from values saved in a file.

    :param name_of_file: Name of the file from which values are loaded.
    :return: 0 if successful
    """
    temperatures, energies, sigmas = load_property_from_file(name_of_file)

    fig, ax = plt.subplots()

    ax.errorbar(temperatures, energies, yerr=sigmas, fmt='ro', label="Measurement")

    ax.set_xlim((0.0, temperatures[-1] + 0.2))
    ax.set_xlabel(r"Temperature $\frac{k_B T}{J}$", fontsize=18)
    ax.set_ylabel(r"Energy per spin $e$", fontsize=18)

    ax.legend(loc="best", fontsize=16)
    ax.grid(visible=True)

    plt.show()
    
    return 0


def save_temp_corr_time(temperature, tau, name_of_file):
    """
    Saves correlation time tau of the system together with the temperature in a file. Saving style: 'temperature tau'.

    :param temperature:(float) Temperature of the system.
    :param tau: (float) Correlation time of the system.
    :param name_of_file: (string) Name of the file in which the values are saved.
    :return: 0 if successful
    """
    with open(name_of_file, "a") as file:
        file.write("%f %f\n" % (temperature, tau))

    return 0

def save_temp_property(temperature, mean, sigma, name_of_file):
    """
    Saves measured mean and error of the given property of the system together with the temperature of the system into a file.
    Saving style: 'temperature mean error'.

    :param temperature: (float) Temperature of the system.
    :param mean: (float) Calculated mean value of the given property of the system (magnetization/energy/susceptibility/specific heat).
    :param sigma: (float) Calculated error of a given property of the system.
    :param name_of_file: (string) Name of the file in which the values are saved.
    :return: 0 if successful
    """
    with open(name_of_file, "a") as file:
        file.write("%f %f %f\n" % (temperature, mean, sigma))

    return 0
