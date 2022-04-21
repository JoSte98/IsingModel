# IsingModel
This repository includes a module to perform a Monte Carlo Simulation of a 2D Grid of Spins interacting via an near neighbor Ising Hamiltonian. 
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
the optional parameter 'init_type' during the creation of the object. To see pre-programed initial states check the docstrings.

Example: Example of creation of a class object and simulation with automatic saving and plotting can be found in the main.py file.

## Simulating statistical properties
To measure statistical properties using minimal number of measurrement blocks N one has to use object.measure_stat_prop(t,N), where t is the number of 
sweeps during which we measure the correlation time of the system. The measured properties are then stored as object.final_magnetization, object.final_energy,
object.final_susceptibility, object.final_specific_heat. If we add _sigma to the end of previously named variables, we get tge stored errors of this properties.
These measured values then can be saved to a file via external plotting functions (see section 'Plotting functions'). 

Note: As for optional parameters, if we want to enable saving of all measured magnetizations and energies to a file, change 'save' to True.

Example: Detailed example of this is the main part of the pre-prepared main.py file. It measures statistical properties for a range of temperatures, 
saves results in the files - NAME_PHASE_FILE (magnetizations), NAME_CORR_FILE (correlation times), NAME_ENERGY_FILE (energies), NAME_SPECIFIC_HEAT_FILE (specific heats),
NAME_SUSCEPTIBILITY_FILE (magnetic susceptibilities) - and then plots them.

## Adding external magnetic field
An external magnetic field h can be added to the simulation during initialization of the IsingModel object via the optional parameter 'external_field'. 
It is by default set to 0.0, but for the sake of clarity we have explicitly initialized it to 0.0 once more in the main.py file to make the example more clear.

## Plotting configurations of the system
A current configuration of the system can be plotted via the function object.plot_state(). 

Note: There is however also an example (not using this built-in function) in the main.py file which includes a piece of code which will produce the state of 
a system after 100 additional sweeps for temperatures in range (1.0,4.0) and plots them.

## Plotting functions
Aside from the in-built plotting functions of the IsingModel class there are also some external functions in plotting.py. Via using 
save_temp_property(temperature, mean, error, name_of_file) we can save any of the previously measured statistical properties in the file (name_of_file: NAME_CORR_FILE,
NAME_PHASE_FILE, NAME_ENERGY_FILE, NAME_SPECIFIC_HEAT_FILE, NAME_CORR_FILE) as the trio of numbers (temperature, mean, error). From this file they then can be loaded 
and plotted for all temperatures (all values for different temperatures saved so far) using functions plot_phase_diagram(NAME_PHASE_FILE), plot_temp_energy(NAME_ENERGY_FILE), 
plot_temp_specific_heat(NAME_SPECIFIC_HEAT_FILE), plot_temp_susceptibility(NAME_SUSCEPTIBILITY_FILE).
Since the correlation time tau does not have an error value, it is due to its different format saved using the function save_temp_corr_time(temperature, tau, name_of_file).
From this file it can then be plotted using plot_temp_corr_time(NAME_CORR_FILE)

Example: In the main.py file the default example contains exactly this saving to file and then plotting of the measured properties for all measured temperatures.
