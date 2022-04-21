"""
Simulation class of a 2D Ising Model via Monte Carlo integration.

@author: Magdalena and Johannes
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.patches as mpatches

class IsingModel:

    def __init__(self, temperature, length=50, init_type='up', equilibrate_steps=200, external_field = 0.0):
        """
        Constructor of the Ising Model class.

        :param temperature: (float) Dimensionless temperature in units of the coupling constant J.
        :param length: (int) Number of spins in the 1D chain. The grid will be of dimension length^2.
        :param init_type: (from {'up', 'down', 'random'}, optional) Type of the initialization of the 2D spin grid, default value 'up'.
        :param equilibrate_steps:(int, optional) Number of sweeps made to equilibrate the system during initialization, default value 200.
        :param external_field:(float, optional) Adds external magnetic field H to the system, defaul value 0.0.
        """
        self.temperature = temperature
        self.length = int(length)
        self.initialize_state(init_type=init_type)
        self.external_field = external_field

        self.magnetizations = []
        self.energies = []
        self.susceptibilities = []
        self.specific_heats = []
        
        self.p_4J_m2H = self.probability(abs(4-2*external_field))
        self.p_8J_m2H = self.probability(abs(8-2*external_field))
        self.p_4J_p2H = self.probability(abs(4+2*external_field))
        self.p_8J_p2H = self.probability(abs(8+2*external_field))
        self.p_p2H = self.probability(abs(2*external_field))
        

        self.equilibrate(equilibrate_steps)

    def initialize_state(self, init_type):
        """
        Initialize a grid of -1 (down) or 1 (up) in the order specified by init_type. Here 'up' means all spins = 1, 'down' means all spins = -1, 'random' random choice of +-1 for each spin. As for 'close_to_opt' means that, depending on the expected magnetization (Onsager formula), a certain number of random spins is switched from 1 to -1.

        :param init_type: (from {'up', 'down', 'random'}) Type of initial state. 
        :return: 0 if successful
        """

        self.state = np.ones((self.length, self.length), dtype=int)

        if init_type == 'up':
            return 0
        elif init_type == 'down':
            self.state *= -1
            return 0
        elif init_type == "close_to_opt":
            # theoretical expected amount of ups to be flipped for expected magnetization of the system
            if self.temperature < 2.369:
                N = int(((1 - (1 - self.temperature/2.369)**(1/8))/2) * self.length**2)
            else:
                N = int(self.length**2/2)
            for i in range(N):
                flip_spin_x = np.random.randint(0, self.length)
                flip_spin_y = np.random.randint(0, self.length)
                self.state[flip_spin_x, flip_spin_y] *= -1
            return 0
        elif init_type == 'random':
            order = np.random.choice([-1, 1], size=int(self.length**2))
            for i in range(self.length):
                for j in range(self.length):
                    self.state[i, j] = order[self.length*i + j]
            return 0
        else: 
            print("This initial state is not known.")

    def equilibrate(self, equilibrate_steps):
        """
        Equilibration of the system for a given number of sweeps, magnetizations and energies are not saved during this procedure.
        
        :param equilibrate_steps: (int) Number of sweeps to equilibrate the system.
        :return: 0 if successful
        """
        for i in range(equilibrate_steps):
            self.simulation_unit(save=False)
        return 0

    def measure_corr_time(self, t_max, plot=False, save_magnetizations=False):
        """
        Evolves the system for t_max sweeps and from attained magnetizations calculates the approximate correlation time of the system.
        
        :param t_max: number of sweeps made in total
        :param plot: (True/False, optional) allows/disallowes plotting of the integrated function chi(t)/chi(0), default value is False
        :param save_magnetizations: (True/False, optional) Allowing/disallowing saving of magnetizations into a file, default value False.
        :return:  0 if successful
        """
        self.simulate(t_max,plot=False, save=save_magnetizations)

        chi = np.zeros(t_max)
        for t in range(t_max):
            m_t_tp = np.array(self.magnetizations[t:t_max])
            m_tp = np.array(self.magnetizations[0: (t_max-t)])

            chi[t] = 1/(t_max - t) * (np.sum(m_tp*m_t_tp) - 1/(t_max - t) * np.sum(m_tp)*np.sum(m_t_tp))

        if plot==True:
            fig, ax = plt.subplots()

            ax.scatter(range(t_max), chi/chi[0], linewidths=1.5, c="r", label="Measurement")

            ax.set_xlabel(r"$t$", fontsize=18)
            ax.set_ylabel(r"$\chi(t)$", fontsize=18)

            ax.legend(loc="best", fontsize=16)
            ax.grid(visible=True)

        self.tau = 0
        for t in range(t_max):
            if chi[t] < 0:
                break
            self.tau += chi[t]
        self.tau /= chi[0]

        return 0

    def measure_stat_prop(self, t_max, num_boxes, plot=False, save_magnetizations=False):
        """
        First measures correlation time of the system, then determines the size of independent measurement blocks. Splits simulated data from correlation time measurement into measurement blocks and calculated corresponding susceptibilities and specific heats. If the number of optained measurements is smaller then num_boxes (minimum required), then additional block are simulated and measured.
        Mean and error of magnetization and energy is then calculated from all measurements done. Susceptibility and specific heat mean and error is calculated from independent block measurements.
        
        :param t_max: (int) Number of sweeps to measure correlation time.
        :param num_boxes: (int) Minimal number of measurement blocks to be performed.
        :param plot: (True/False, optional) Allows/disallows plotting during measuring of the correlation time, default value False.
        :param save_magnetizations: (True/False, optional) Allowing/disallowing saving of magnetizations into a file, default value False.
        :return: 0 if successful
        """
        self.measure_corr_time(t_max, plot, save_magnetizations)

        tau = min(200, self.tau)
        t_box = int(16*tau)
        num_boxes_done = math.floor(t_max/t_box) #how many blocks do we have already
        num_boxes_needed = num_boxes - num_boxes_done #how many blocks need to be additionaly measured
        N = int(max(0, t_box - t_max)) #how much steps to round the next block

        for i in range(num_boxes_done):
            start = int((i)*t_box)
            stop = int((i +1)*t_box)
            self.susceptibilities.append(self.susceptibility(self.magnetizations[start:stop]))
            self.specific_heats.append(self.specific_heat(self.energies[start:stop]))

        if num_boxes_needed>0:
            self.simulate(N, plot=False,save=save_magnetizations)
            for i in range(num_boxes_needed):
                self.simulate(t_box, plot=False,save=save_magnetizations)
                self.susceptibilities.append(self.susceptibility(self.magnetizations[-t_box:]))
                self.specific_heats.append(self.specific_heat(self.energies[-t_box:]))

        #means calculation
        self.final_susceptibility = np.mean(self.susceptibilities)
        self.final_specific_heat = np.mean(self.specific_heats)
        self.final_magnetization = np.mean(np.absolute(self.magnetizations))
        self.final_energy = np.mean(self.energies)/self.length**2
        
        #errors calculation
        self.final_susceptibility_sigma = np.sqrt( ( np.mean(np.array(self.susceptibilities)**2) -
                                                     np.mean(np.array(self.susceptibilities))**2 ) )
        self.final_specific_heat_sigma = np.sqrt( ( np.mean(np.array(self.specific_heats)**2) -
                                                     np.mean(np.array(self.specific_heats))**2 ) )
        self.final_magnetization_sigma = np.sqrt( 2*tau/(t_max + N + (num_boxes_needed)*16*tau) *
                                                   ( np.mean(np.absolute(self.magnetizations)**2) -
                                                     np.mean(np.absolute(self.magnetizations))**2 ) )
        self.final_energy_sigma = np.sqrt( 2*tau/(t_max + N + (num_boxes_needed)*16*tau) *
                                                   ( np.mean(np.array(self.energies)**2) -
                                                     np.mean(np.array(self.energies))**2 ) )/self.length

        return 0

    def susceptibility(self, magnetizations):
        """
        Calculates susceptibility from a given list of measured magnetizations.

        :param magnetizations: (list of floats) List of average magnetizations per spin from a given measurement block.
        :return: Susceptibility.
        """
        susceptibility = 1/self.temperature * (np.mean(np.array(magnetizations)**2) - np.mean(np.array(magnetizations))**2) *self.length**2
        return susceptibility

    def specific_heat(self, energies):
        """
        Calculates specific heat from a given list of measured total energies.

        :param E: (list of floats) List of total energies of the system from a given measurement block
        :return: Specific heat as a float.
        """
        C = 1/(self.length**2 * self.temperature**2) * (np.mean(np.array(energies)**2) - np.mean(np.array(energies))**2)
        return C

    def new_state(self):
        """
        Creates a new state out of the current state via flipping randomly chosen spin and half of the energy difference of the system caused by the flip.

        :return: New state, half of the corresponding change in total energy of the system.
        """
        flip_spin_x = np.random.randint(0, self.length)
        flip_spin_y = np.random.randint(0, self.length)

        new_state = self.state.copy()
        
        dE = (self.state[flip_spin_x, (flip_spin_y+1)%self.length] + self.state[flip_spin_x, (flip_spin_y-1)%self.length]\
             + self.state[(flip_spin_x+1)%self.length, flip_spin_y] + self.state[(flip_spin_x-1)%self.length, flip_spin_y])\
             * self.state[flip_spin_x, flip_spin_y] + self.external_field*new_state[flip_spin_x, flip_spin_y]

        new_state[flip_spin_x, flip_spin_y] *= -1
        

        return new_state, dE

    def energy(self, state):
        """
        Calculates the energy (Ising hamiltonian) of the state in dimensionless units J/(k_B T).

        :param state: (state of the system - 2D grid) State, whose energy should be calculated.
        :return: Energy as a float.
        """

        energy = 0.0
        for i in range(self.length):
            for j in range(self.length):
                # left and right neighbor
                energy -= state[i, j]*(state[i, (j+1)%self.length] + state[i, (j-1)%self.length])
                # top and bottom neighbor
                energy -= state[i, j]*(state[(i+1)%self.length, j] + state[(i-1)%self.length, j])
                energy -=  2* self.external_field*state[i,j]
                

        energy *= 0.5

        return energy

    def probability(self, energy):
        """
        Calculates the probability of accepting new state of the Ising system for a given total energy difference caused by this change.

        :param energy: (float) Energy difference (energy of the current state of the system - energy of the proposed state).
        :return: Probability as float.
        """

        return np.exp(-energy/self.temperature)

    def magnetization(self):
        """
        Measures the average spin (magnetization) of the current state of the system.

        :return: Magnetization as float.
        """

        M = 1/(self.length**2) * np.sum(self.state)

        return M

    def evolve(self):
        """
        Draws new state by flipping randomly chosen spin and depending on the change of the total energy of the system accepts the new state with a given probability.
        
        :return: 0 if successful.
        """
        new_state , dE = self.new_state()

        if dE<=0: #energy is lowering/staying the same -> automaticaly accept new state
            self.state = new_state

        else: #accepting the state depending on the energy difference if the new energy is bigger
            if (abs(dE-abs(4.00-self.external_field))<10e-8):
                p = self.p_8J_m2H
                #print(p)
            elif (abs(dE-abs(2.00-self.external_field))<10e-8):
                p = self.p_4J_m2H
            elif (abs(dE-abs(4.00+self.external_field))<10e-8):
                p = self.p_8J_p2H
            elif (abs(dE-abs(2.00+self.external_field))<10e-8):
                p = self.p_4J_p2H
            else:
                p = self.p_p2H
                
            t = np.random.random()
            if (t<p): # accept the new state
                self.state = new_state

        return 0

    def simulation_unit(self, save=True):
        """
        Does one sweep through the lattice, i. e. tries to swap a random spin and evolve the system self.lenght**2 times (to give every spin a chance to flip).

        :param save: (True/False, optional) Gives an option to turn off (False) measuring of energy and magnetization after every sweep, default value True.   
        :return: 0 if successful
        """
        for i in range(self.length**2):
            self.evolve()

        if save:
            self.energies.append(self.energy(self.state))
            self.magnetizations.append(self.magnetization())

        return 0

    def simulate(self, num_repetitions, plot=True, save=True):
        """
        Does num_repetitions of sweeps through the lattice and saves the magnetizations into a file.
        
        :param num_repetitions: Number of sweeps to be done.
        :param plot: (True/False, optional) Allowing/disallowing the plotting of magnetizations and energies, default value True.
        :param save: (True/False, optional) Allowing/disallowing saving of magnetizations into a file, default value True.
        :return: 0 if successful
        """
        for i in range(num_repetitions):
            self.simulation_unit()
        
        if (save==True):
            self.save_magnetization()
        if (plot == True):
            self.plot_magnetization()
            self.plot_energies()

        return 0

    def plot_magnetization(self):
        """
        Plots average magnetizations per spin of the system saved so far.

        :return: 0 if successful
        """
        fig, ax = plt.subplots()

        x_range = range(len(self.magnetizations))
        ax.scatter(x_range, self.magnetizations, linewidths=2.0)

        ax.set_xlabel("Measurement", fontsize=18)
        ax.set_ylabel("Average Magnetization", fontsize=18)

        ax.set_xticks(x_range)

        ax.set_title("Result Measurements Markov Chain", fontsize=20)

        plt.show()
        
        return 0

    def save_magnetization(self):
        """
        Saves magnetizations measured so far into an txt file of automatically generated name.

        :return: 0 if successful
        """
        num_repetitions = len(self.magnetizations)
        name_of_file = "magn_" + str(self.temperature).replace(".", "") + "_size_" + str(int(self.length)) + ".txt"
        with open(name_of_file, "a") as file:
            file.write("%f %d %d\n" % (self.temperature, num_repetitions, self.length))
            for i in range(num_repetitions):
                file.write("%f " % self.magnetizations[i])
                file.write("\n")

        return 0

    def plot_energies(self):
        """
        Plots total energies of the system measured so far.

        :return: 0 if successful
        """
        fig, ax = plt.subplots()

        x_range = range(len(self.energies))
        ax.plot(x_range, self.energies, linewidth=1.0)

        ax.set_xlabel("Measurement", fontsize=18)
        ax.set_ylabel("Energy", fontsize=18)

        ax.set_xticks(x_range)

        ax.set_title("Result Measurements Markov Chain", fontsize=20)

        plt.show()
        
        return 0

    def plot_state(self):
        """
        Plots the current state of the system (values of spins in the lattice) using imshow.

        :return: 0 if successful
        """
        fig, ax = plt.subplots()
        im = ax.imshow(self.state)
        ax.set_title('Current configuration of the system')
        colors = [ im.cmap(im.norm(1)),im.cmap(im.norm(-1))]
        # create a patch (proxy artist) for every color 
        patches = [ mpatches.Patch(color=colors[0], label="up"),mpatches.Patch(color=colors[1], label="down")]
        # put those patched as legend-handles into the legend
        plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )

        plt.show()
        
        return 0
        
        








