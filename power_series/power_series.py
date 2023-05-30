import numpy as np
import matplotlib.pyplot as plt

"""
a_0 = eval(input('a_0 = '))
a_1 = eval(input('a_1 = '))
epsilon = eval(input('epsilon = '))
shift_xi = eval(input('shift_xi = '))
"""
# make a function that takes previous value and outputs new term
# make another function to call the coeff.


def recursion_relationship(recursion_index, recursion_coefficients, shifted_energy, potential_shift):
    """
    :param recursion_index: n, integer
    :param recursion_coefficients: a_n-2, a_n-1, a_n, list of floats
    :param shifted_energy: λ, float
    :param potential_shift: ±ξ0, float
    :return new_coefficient: a_n+2
    """
    new_coefficient = (- 2 * shifted_energy * recursion_coefficients[2] + recursion_coefficients[0] +
                         2 * potential_shift * recursion_coefficients[1]) /\
                      ((recursion_index + 1)*(recursion_index + 2))
    return new_coefficient


def wave_function(reduced_position, maximum_coefficient, starting_coefficients, simulation_parameters):
    """
    :param reduced_position: ξ
    :param maximum_coefficient: n_max
    :param starting_coefficients: a0, a1, a2. a3
    :param simulation_parameters: λ, ξ0
    :return: ψ(ξ)
    """
    coefficients = starting_coefficients
    for index in range(len(starting_coefficients), maximum_coefficient+1):
        coefficients.append(recursion_relationship(index, coefficients[index-4:index-1],
                                                   simulation_parameters[0], simulation_parameters[1]))
    return np.polynomial.polynomial.polyval(reduced_position, coefficients)


def indicial_equations(eigenvalue):
    a_2 = -eigenvalue*a_0
    if shift_xi >= 0:
        a_3 = -(eigenvalue*a_1 - shift_xi*a_0)/3
    if shift_xi < 0:
        a_3 = -(eigenvalue*a_1 + shift_xi*a_0)/3

    #if shift_xi == 0:
    #    a_3 = 'this is a harmonic oscillator'
    return [a_2, a_3]


def recursion_relation(indicial_terms, eigenvalue):
    a_2, a_3 = indicial_terms
    if shift_xi >= 0:
        a_4 = (-2*eigenvalue*a_2 + a_0 + 2*shift_xi*a_1)/(4*3)
        a_5 = (-2*eigenvalue*a_3 + a_1 + 2*shift_xi*a_2)/(5*4)
        a_6 = (-2*eigenvalue*a_4 + a_2 + 2*shift_xi*a_3)/(6*5)
        a_7 = (-2*eigenvalue*a_5 + a_3 + 2*shift_xi*a_4)/(7*6)
        a_8 = (-2*eigenvalue*a_6 + a_4 + 2*shift_xi*a_5)/(8*7)
        a_9 = (-2*eigenvalue*a_7 + a_5 + 2*shift_xi*a_6)/(9*8)
        a_10 = (-2*eigenvalue*a_8 + a_6 + 2*shift_xi*a_7)/(10*9)
        # a_terms = [a_0, a_1, indicial_terms[0], indicial_terms[1], a_4, a_5, a_6, a_7, a_8, a_9]
    if shift_xi < 0:
        a_4 = (-2*eigenvalue * a_2 + a_0 - 2 * shift_xi * a_1) / (4 * 3)
        a_5 = (-2*eigenvalue * a_3 + a_1 - 2 * shift_xi * a_2) / (5 * 4)
        a_6 = (-2*eigenvalue * a_4 + a_2 - 2 * shift_xi * a_3) / (6 * 5)
        a_7 = (-2*eigenvalue * a_5 + a_3 - 2 * shift_xi * a_4) / (7 * 6)
        a_8 = (-2*eigenvalue * a_6 + a_4 - 2 * shift_xi * a_5) / (8 * 7)
        a_9 = (-2*eigenvalue * a_7 + a_5 - 2 * shift_xi * a_6) / (9 * 8)
        a_10 = (-2*eigenvalue*a_8 + a_6 - 2*shift_xi*a_7)/(10*9)

    a_terms = [a_0, a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8, a_9, a_10]
    return a_terms


def plot_expansion(coefficients):
    xi = np.linspace(-10, 10, 1000)
    psi_array = []
    # for i coef do a_i*xi**i
    for i in xi:
        psi = coefficients[0] * i**0 + coefficients[1] * i**1 + coefficients[2] * i**2 + coefficients[3] * i**3 + \
              coefficients[4] * i**4 + coefficients[5] * i**5 + coefficients[6] * i**6 + coefficients[7] * i**7 + \
              coefficients[8] * i**8 + coefficients[9] * i**9 + coefficients[10] * i**10
        psi_array.append(psi)
    x_label = plt.xlabel('xi')
    y_label = plt.ylabel('psi')
    plt.plot(xi, psi_array)
    plt.ylim([-1, 1])
    plt.show()
    return psi_array


if __name__ == '__main__':
    """
    a_0 = 1, a_1 = 0, shift = 0, eigenvalue = 0.5
    a_0 = 0, a_1 = 1, shift = 0, eigenvalue = 1.5
    a_0 = 1, a_1 = 0, shift = 0, eigenvalue = 2.5
    .
    .
    .
    """
    a_0 = 0
    a_1 = 1
    # eigenvalue should be 0.5, (bell curve)
    # then set a_0 = 0, a_1 = 1, and search for another value
    eigenvalue = 1.5   # (epsilon - shift_xi^2/2)
    shift_xi = 0
    epsilon = eigenvalue + shift_xi**2/2
    print(f'a_2 = {indicial_equations(eigenvalue)[0]},'
          f'a_3 = {indicial_equations(eigenvalue)[1]}')
    # print(f' a_terms = {recursion_relation(indicial_equations(eigenvalue), eigenvalue)}')
    #print(plot_expansion(recursion_relation(indicial_equations(eigenvalue), eigenvalue)))
    initial_coefficients = [a_0, a_1]
    for coefficient in indicial_equations(eigenvalue):
        initial_coefficients.append(coefficient)

    reduced_positions = np.linspace(-10, 10, num=1_000)
    maximum_index = 99
    psi_values = wave_function(reduced_positions, maximum_index, initial_coefficients, [eigenvalue, shift_xi])

    plt.plot(reduced_positions, psi_values)
    # plt.xlim([-10, 10])
    # plt.ylim([-1, 1.5])
    plt.show()

    # print(epsilon)
