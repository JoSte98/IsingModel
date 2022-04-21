# IsingModel
This repository includes a module to simulate Monte Carlo Simulation of a 2D Grid of Spins interacting via an near neighbor Ising Hamiltonian. 
It allows to compute the statistical properties of the system such as magnetization, averagre energy per spin, magnetic susceptibility and 
specific heat of the system. 
By simply pulling the repository you should be able to execute the main.py file. This file contains different important parts of the
module: First the main class “IsingModel” is importet via “from IsingModel import IsingModel” and also external plotting functions via 
"from plotting import *". 
Executing the main.py file in it's default form you will simulate the system for tempreratures in ranges (1.0,4.0), save and plot all
statictical properties and correlation time graphs. Here we give some basic overwiev of a few most important functions and also explain
what other examples can be found in the main.py file.

## Simulation
The basis of every simulation would be the creation of an object of class IsingModel with specified temperature and the size of the lattice. 
After that, one can run a simulation of the system for a specified number of sweeps t via object.simulate(t). This function would also automaticaly
save to a file and plot the magnetizations and energies of the system collected during the simulation, these can be disablet via setting optional 
parameters save=False and plot=False. These properties then can also be plotted additionaly after simulation via object.plot_magnetizations() and 
object.plot_energies().

Note: The system automatically initialize in the configuration 'up' meaning all spins are +1, there is a few other options which can be accessed via 
the optional parameter 'init_type' during the creation of the object. To see pre-programed initial states check docstrings.

Example: Example of creation of an class object and simulation with automatic saving and plotting can be found in the main.py file.

## Simulating statistical properties


Example: Detailed example of this is the main part of the pre-prepared main.py file. It measures statistical properties for a range of temperatures, 
saves results in the files - NAME_PHASE_FILE (magnetizations), NAME_CORR_FILE(correlation times), NAME_ENERGY_FILE (energies), NAME_SPECIFIC_HEAT_FILE(specific heats),
NAME_SUSCEPTIBILITY_FILE (magnetic susceptibilities) - and then plots them.

## Adding external magnetic field


## Plotting configurations of the system



## Plotting functions