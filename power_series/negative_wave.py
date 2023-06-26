import numpy as np
import matplotlib.pyplot as plt


def negative_wave_function(starting_coefficients, xi_array):
    """
    :param starting_coefficients: [a0, a1, a2, a3] list of floats
    :param xi_array: ξ array of values
    :return:
    """
    maximum_index = 10_000  # y values kinda
    coefficients = starting_coefficients
    for index in range(len(starting_coefficients), maximum_index + 1):
        coefficients.append(
            recurence_relation_neg(coefficients[index - 4:index - 1], shifted_xi, shifted_eigenvalue, index)
        )
    return np.polynomial.polynomial.polyval(xi_array, coefficients)


def recurence_relation_neg(starting_coefficients, shifted_xi, shifted_eigenvalue, recurence_index):
    """
    :param starting_coefficients: [a_n-2, a_n-1, a_n] list of floats
    :param shifted_xi: ξ0, float
    :param shifted_eigenvalue: λ, float
    :param recurence_index: starting at n = 2, thus getting a4, integer
    :return: new_coefficient: a_n+2
    """
    new_coefficient = starting_coefficients[0] - 2 * shifted_xi * starting_coefficients[1] \
                      - 2 * starting_coefficients[2] * shifted_eigenvalue / \
                      (recurence_index * (recurence_index - 1))
    return new_coefficient


def indicial_equations_neg(a0, a1, shifted_xi, shifted_eigenvalue):
    """
    :param a0: float
    :param a1: float
    :param shifted_xi: float
    :param shifted_eigenvalue: float
    :return: 2 floats in a list
    """
    a2 = -shifted_eigenvalue*a0
    a3 = -1/3*shifted_eigenvalue*a1 - 1/3*shifted_xi*a0
    return [a0, a1, a2, a3]


if __name__ == "__main__":
    a0 = 1
    a1 = 0
    shifted_xi = 0
    epsilon = 0.5
    shifted_eigenvalue = epsilon - shifted_xi**2/2

    #print(type(indicial_equations_neg(a0, a1, shifted_xi, shifted_eigenvalue)))
    starting_coefficients = indicial_equations_neg(a0, a1, shifted_xi, shifted_eigenvalue)
    #print(type(starting_coefficients))

    number_of_positions = 10_001  # x values kinda
    xi_array = np.linspace(-10, 0, num=number_of_positions)
    ##
    negative_psi_values = negative_wave_function(starting_coefficients, xi_array)
    plt.plot(xi_array, negative_psi_values)
    plt.xlim([-10, 10])
    plt.ylim([-3.5, 3.5])
    plt.axhline()
    plt.axvline(0)
    plt.show()
