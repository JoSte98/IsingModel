"""
Simulation class of a 2D Ising Model via Monte Carlo integration.

@author: Magdalena and Johannes
"""

import numpy as np
import matplotlib.pyplot as plt

class IsingModel:

    def __init__(self, temperature, length=50, init_type='up', equilibrate_steps=200):
        """
        Constructor of the Ising Model class.

        :param temperature: Dimensionless temperature in units of the coupling constant J.
        :param length: Number of spins in the 1D chain. The grid will be of dimension length^2.
        :param init_type: From {'up', 'down', 'random'}. Type of the initialization of the 2D spin grid.
        :param equilibrate_steps: Number of weeps to equilibrate the system.
        """
        self.temperature = temperature
        self.length = int(length)
        self.initialize_state(init_type=init_type)

        self.magnetizations = []
        self.energies = []
        self.susceptibilities = []
        self.specific_heats = []
        
        self.p_4J = self.probability(4)
        self.p_8J = self.probability(8)

        self.equilibrate(equilibrate_steps)

    def initialize_state(self, init_type):
        """
        Initialize a grid of -1 (down) or 1 (up) in the order specified by init_type.

        :param init_type: From {'up', 'down', 'random'}, where up means all spins = 1, down = -1, and random random.
        :return: 0 if successfull
        """

        self.state = np.ones((self.length, self.length), dtype=int)

        if init_type == 'up':
            return 0
        elif init_type == 'down':
            self.state *= -1
            return 0
        elif init_type == "close_to_opt":
            # theoretical expected amount of ups to be flipped for perfect magnetization
            if self.temperature < 2.369:
                N = int(((1 - (1 - self.temperature/2.369)**(1/8))/2) * self.length**2)
            else:
                N = int(self.length**2/2)
            print(N/self.length**2)
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
            #self.plot_state()
            return 0
        elif init_type == 'half-up-half-down':
            for i in range(self.length):
                for j in range(self.length/2):
                    self.state[i,j] *= -1
            return 0
        else: 
            print("This initial state is not known.")

    def equilibrate(self, equilibrate_steps):
        """

        :param equilibrate_steps:
        :return: 0 if successful
        """
        for i in range(equilibrate_steps):
            self.simulation_unit(save=False)
        return 0

    def measure_corr_time(self, t_max, plot=False):
        """

        :param t_max:
        :return:
        """
        for t in range(t_max):
            self.simulation_unit()

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

    def measure_stat_prop(self, t_max, num_boxes, plot=False):
        """
        :return:
        """
        self.measure_corr_time(t_max, plot)

        tau = min(100, self.tau)

        N = max(0, 16*tau - t_max)
        N = int(N)
        self.simulate(N, plot=False)
        self.susceptibilities.append(self.susceptibility(self.magnetizations))
        self.specific_heats.append(self.specific_heat(self.energies))

        for i in range(num_boxes - 1):
            self.simulate(int(16*tau), plot=False)
            self.susceptibilities.append(self.susceptibility(self.magnetizations[-int(16*tau):]))
            self.specific_heats.append(self.specific_heat(self.energies[-int(16*tau):]))

        self.final_susceptibility = np.mean(self.susceptibilities)
        self.final_specific_heat = np.mean(self.specific_heats)
        self.final_magnetization = np.mean(self.magnetizations)
        self.final_energy = np.mean(self.energies)

        self.final_susceptibility_sigma = np.sqrt( 2*tau/(num_boxes*t_max) *
                                                   ( np.mean(np.array(self.susceptibilities)**2) -
                                                     np.mean(np.array(self.susceptibilities))**2 ) )
        self.final_specific_heat_sigma = np.sqrt( 2*tau/(num_boxes*t_max) *
                                                   ( np.mean(np.array(self.specific_heats)**2) -
                                                     np.mean(np.array(self.specific_heats))**2 ) )
        self.final_magnetization_sigma = np.sqrt( 2*tau/(num_boxes*t_max) *
                                                   ( np.mean(np.array(self.magnetizations)**2) -
                                                     np.mean(np.array(self.magnetizations))**2 ) )
        self.final_energy_sigma = np.sqrt( 2*tau/(num_boxes*t_max) *
                                                   ( np.mean(np.array(self.energies)**2) -
                                                     np.mean(np.array(self.energies))**2 ) )

        return 0

    def susceptibility(self, m):
        """

        :param m:
        :return:
        """

        return 1/self.temperature * (np.sum(np.array(m)**2) - np.sum(np.array(m))**2)

    def specific_heat(self, E):
        """

        :param E:
        :return:
        """
        return 1/(self.length**2 * self.temperature**2) * (np.sum(np.array(E)**2) - np.sum(np.array(E))**2)

    def new_state(self):
        """
        Creates out of the current state a new one via flipping a random spin.

        :return: New state.
        """
        flip_spin_x = np.random.randint(0, self.length)
        flip_spin_y = np.random.randint(0, self.length)

        new_state = self.state.copy()
        
        dE = (self.state[flip_spin_x, (flip_spin_y+1)%self.length] + self.state[flip_spin_x, (flip_spin_y-1)%self.length]\
             + self.state[(flip_spin_x+1)%self.length, flip_spin_y] + self.state[(flip_spin_x-1)%self.length, flip_spin_y])\
             * self.state[flip_spin_x, flip_spin_y]
        #print(dE)
        new_state[flip_spin_x, flip_spin_y] *= -1
        

        return new_state, dE

    def energy(self, state):
        """
        Calculates the energy (Ising hamiltonian) of the state in dimensionless units J/(k_B T).

        :param state: State, whose energy should be calculated (2D grid).
        :return: Energy as a float.
        """

        energy = 0.0
        for i in range(self.length):
            for j in range(self.length):
                # left and right neighbor
                energy -= state[i, j]*(state[i, (j+1)%self.length] + state[i, (j-1)%self.length])
                # top and bottom neighbor
                energy -= state[i, j]*(state[(i+1)%self.length, j] + state[(i-1)%self.length, j])

        energy *= 0.5

        return energy

    def probability(self, energy):
        """
        Calculates the probability of a state in the Ising system.

        :param state: State, whose energy should be calculated (2D grid).
        :return: Probability as float.
        """

        return np.exp(-energy/self.temperature)

    def magnetization(self):
        """
        Measures the normed magnetization of the current state.

        :return: Magnetization as float.
        """

        M = 1/(self.length**2) * np.sum(self.state)

        return M

    def evolve(self):
        """

        :return:
        """
        new_state , dE = self.new_state()
        #new_state_energy = self.energies[-1]+2*dE

        if dE<=0:
            #print("right away")
            self.state = new_state
            #self.energies.append(new_state_energy)
            #print("tried")
        else:
            if (abs(dE-4.00)<10e-8):
                p = self.p_8J
                #print(p)
            else:
                p = self.p_4J
                #print(p)
            t = np.random.random()
            if (t<p): # accept the new state
                self.state = new_state
                #self.energies.append(new_state_energy)
                #print("accepted")

            

        return 0

    def simulation_unit(self, save=True):
        """

        :return:
        """
        for i in range(self.length**2):
            self.evolve()

        if save:
            self.energies.append(self.energy(self.state))
            self.magnetizations.append(self.magnetization())

        return 0

    def simulate(self, num_repetitions, plot=True):
        """

        :param num_repetitions:
        :return:
        """
        for i in range(num_repetitions):
            self.simulation_unit()

        self.save_magnetization()
        if (plot == True):
            self.plot_magnetization()

        return 0

    def plot_magnetization(self):
        fig, ax = plt.subplots()

        x_range = range(len(self.magnetizations))
        ax.scatter(x_range, self.magnetizations, linewidths=2.0)

        ax.set_xlabel("Measurement", fontsize=18)
        ax.set_ylabel("Average Magnetization", fontsize=18)

        ax.set_xticks(x_range)

        ax.set_title("Result Measurements Markov Chain", fontsize=20)

        plt.show()

    def save_magnetization(self):
        """

        :return:
        """
        num_repetitions = len(self.magnetizations)
        name_of_file = "magn_" + str(self.temperature).replace(".", "") + "_size_" + str(int(self.length)) + ".txt"
        with open(name_of_file, "w") as file:
            file.write("%f %d %d\n" % (self.temperature, num_repetitions, self.length))
            for i in range(num_repetitions):
                file.write("%f " % self.magnetizations[i])
                file.write("\n")

        return 0

    def plot_energies(self):
        """

        :return:
        """
        fig, ax = plt.subplots()

        x_range = range(len(self.energies))
        ax.plot(x_range, self.energies, linewidth=1.0)

        ax.set_xlabel("Measurement", fontsize=18)
        ax.set_ylabel("Energy", fontsize=18)

        ax.set_xticks(x_range)

        ax.set_title("Result Measurements Markov Chain", fontsize=20)

        plt.show()

    def plot_state(self):
        """

        :return:
        """
        fig, ax = plt.subplots()
        ax.imshow(self.state)
        plt.show()
        return 0
        
        








