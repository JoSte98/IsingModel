"""
Simulation class of a 2D Ising Model via Monte Carlo integration.

@author: Magdalena and Johannes
"""

class IsingModel:

    def __init__(self, temperature, length=50, init_type='up'):
        """
        Constructor of the Ising Model class.

        :param temperature: Dimensionless temperature in units of the coupling constant J.
        :param length: Number of spins in the 1D chain. The grid will be of dimension length^2.
        :param init_type: From {'up', 'down', 'random'}. Type of the initialization of the 2D spin grid.
        """
        self.temperature = temperature
        self.length = int(length)
        self.initialize_state(init_type=init_type)

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
        elif init_type == 'random':
            order = np.random.choices([-1, 1], k=int(self.length**2))
            for i in range(self.length):
                for j in range(self.length):
                    self.state[i, j] = order[self.length*i + j]
            return 0

    def new_state(self):
        """
        Creates out of the current state a new one via flipping a random spin.

        :return: New state.
        """
        flip_spin_x = np.random.randint(0, self.length-1)
        flip_spin_y = np.random.randint(0, self.length-1)

        new_state = self.state.copy()

        new_state[flip_spin_x, flip_spin_y] *= -1

        return new_state

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

        energy *= 0.5/self.temperature

        return energy

    def probability(self, state):
        """
        Calculates the probability of a state in the Ising system.

        :param state: State, whose energy should be calculated (2D grid).
        :return: Probability as float.
        """

        return np.exp(-self.energy(state))

    def magnetization(self):
        """
        Measures the normed magnetization of the current state.

        :return: Magnetization as float.
        """

        M = 1/(self.length**2) * np.sum(self.state)

        return M






