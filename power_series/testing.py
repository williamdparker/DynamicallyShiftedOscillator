import numpy as np
import matplotlib.pyplot as plt


# input 4 [a0 a1 a2 a3] lambda shift_xi
def recursion_relationship(recursion_index, recursion_coefficients, shifted_energy, potential_shift):
    """
    :param recursion_index: n+2, integer
    :param recursion_coefficients: a_n-2, a_n-1, a_n, list of floats
    :param shifted_energy: λ, float
    :param potential_shift: ±ξ0, float
    :return new_coefficient: a_n+2
    """
    new_coefficient = (- 2 * shifted_energy * recursion_coefficients[2] + recursion_coefficients[0] +
                       2 * potential_shift * recursion_coefficients[1]) / \
                      ((recursion_index - 1) * (recursion_index))

    # new_coefficient = -recursion_coefficients[2]/(recursion_index+2)
    return new_coefficient  # a4


def wave_function(reduced_position, maximum_coefficient, starting_coefficients, simulation_parameters):
    """
    :param reduced_position: ξ
    :param maximum_coefficient: n_max
    :param starting_coefficients: a0, a1, a2. a3
    :param simulation_parameters: λ, ξ0
    :return: ψ(ξ)
    """
    coefficients = starting_coefficients
    for index in range(len(starting_coefficients), maximum_coefficient + 1):
        coefficients.append(recursion_relationship(index, coefficients[index - 4:index - 1],
                                                   simulation_parameters[0], simulation_parameters[1]))
    return np.polynomial.polynomial.polyval(reduced_position, coefficients)


def indicial_equations(eigenvalue, coefficients):
    a_2 = -eigenvalue * coefficients[0]
    if shift_xi >= 0:
        a_3 = -(eigenvalue * coefficients[1] - shift_xi * coefficients[0]) / 3
    if shift_xi < 0:
        a_3 = -(eigenvalue * coefficients[1] + shift_xi * coefficients[0]) / 3

    # if shift_xi == 0:
    #    a_3 = 'this is a harmonic oscillator'
    return [a_2, a_3]


if __name__ == '__main__':
    a_0 = 1
    a_1 = 0
    shift_xi = 1
    # parkers epsilon: 0.3094688	0.73422037	1.49997346	2.19740304	2.99866356

    #epsilon = 1.6  # 1.500 n=1 shift=1 2.998601 # 6.27841
    #for epsilon in np.linspace(0, 1.5, 15):
    for epsilon in np.linspace(4.3, 4.4, 10):
        # 0.309459 first state (bell curve blown up at the end) note due to no trunkation, positive xi blows up exponentially

        # 1.500000 second state (end behavior again messed up dueto no trunkation)

        eigenvalue = epsilon - shift_xi ** 2 / 2
        """
        Here we wish to loop over many epsilons to find multiply eigenvalues via wag the dog

        for epsilon in range(0, 4, 1):
            eigenvalue = epsilon - shift_xi ** 2 / 2

            reduced_positions = np.linspace(-10, 10, num=1_000)
            maximum_index = 10_000
            initial_coefficients = [a_0, a_1]

            psi_values = wave_function(reduced_positions, maximum_index, initial_coefficients, [eigenvalue, shift_xi])
            recursion_relationship(4, initial_coefficients, eigenvalue, shift_xi)
            plt.plot(reduced_positions, psi_values)
            plt.xlim([-10, 10])
            plt.ylim([-3.5, 3.5])
        plt.show()
    """

        #print(f'a_2 = {indicial_equations(eigenvalue, [a_0, a_1])[0]},'
        #      f'a_3 = {indicial_equations(eigenvalue, [a_0, a_1])[1]}')
        # print(f' a_terms = {recursion_relation(indicial_equations(eigenvalue), eigenvalue)}')
        # print(plot_expansion(recursion_relation(indicial_equations(eigenvalue), eigenvalue)))
        initial_coefficients = [a_0, a_1]
        for coefficient in indicial_equations(eigenvalue, initial_coefficients):
            initial_coefficients.append(coefficient)
        #print(initial_coefficients)
        #print(recursion_relationship(4, initial_coefficients, eigenvalue, shift_xi))  # I agree

        reduced_positions = np.linspace(-10, 10, num=10_001)
        maximum_index = 10_000
        #
        # print(len(reduced_positions)) = 1000
        # print(len(wave_function(reduced_positions, maximum_index, initial_coefficients, [eigenvalue, shift_xi])[1])) = 10,001

        psi_values = wave_function(reduced_positions, maximum_index, initial_coefficients, [eigenvalue, shift_xi])
        #
        plt.plot(reduced_positions, psi_values)
        plt.xlim([-10, 10])
        plt.ylim([-3.5, 3.5])
        plt.axhline()
        plt.axvline(1)
        # plt.savefig("eigenvalue_1")
    plt.show()
        # print(epsilon)

