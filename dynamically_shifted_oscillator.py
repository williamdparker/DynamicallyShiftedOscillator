import numpy as np
import matplotlib.pyplot as plt


def draw_dynamically_shifted_oscillator_potential(function_positions, dynamical_shift, potential_width=1.0,
                                       maximum_potential_value=2.0, plot_axes=np.array([])):
    """
    :param plot_axes:
    :param function_positions:  x-values in a Numpy array
    :param harmonic_oscillator_width:  a
    :param maximum_potential_value: sets maximum V (and minimum V = -maximum V)
    :return:
    """
    force_constant = 1. / potential_width
    minimum_position = np.min(function_positions)
    minimum_potential_position = minimum_position - 0.333 * potential_width
    maximum_position = np.max(function_positions)
    maximum_potential_position = maximum_position + 0.333 * potential_width

    potential_positions = np.linspace(minimum_potential_position, maximum_potential_position, num=10000)
    potential_values = np.piecewise(potential_positions, [potential_positions < 0, potential_positions >= 0],
                        [lambda x: 0.5*force_constant*(x + dynamical_shift)**2,
                         lambda x: 0.5*force_constant*(x - dynamical_shift)**2])

    if plot_axes.size != 0:
        number_of_subplots = np.sum(plot_axes.shape)
        for axis_index in np.arange(number_of_subplots):
            plot_axes[axis_index].fill_between(potential_positions, potential_values,
                                               facecolor='gray', alpha=0.2)
            plot_axes[axis_index].fill_between(potential_positions,
                                               np.full(len(potential_positions), -maximum_potential_value),
                                               facecolor='gray', alpha=0.2)
            plot_axes[axis_index].set_xlim([np.min(potential_positions), np.max(potential_positions)])
            plot_axes[axis_index].set_ylim([-maximum_potential_value, maximum_potential_value])

    else:
        plt.fill_between(potential_positions, potential_values,
                         facecolor='gray', alpha=0.2)
        plt.fill_between(potential_positions, np.full(len(potential_positions), -maximum_potential_value),
                         facecolor='gray', alpha=0.2)

        plt.xlim([np.min(potential_positions), np.max(potential_positions)])
        plt.ylim([-0.5*maximum_potential_value, 0.5*maximum_potential_value])

    return


def dynamically_shifted_oscillator_differential_equation(xi, psi, epsilon, xi_naught):
    """
    Sets up second-order differential equation for solve_ivp as two linear differential equations
    psi' = psi'
    psi'' = ((xi±xi0)**2 - 2*epsilon)*psi
    :param xi:  position value
    :param psi:  list of wave function and wave-function derivative's values
    :param epsilon: adjustable reduced energy parameter to find solutions
    :param xi_naught: dynamical shift
    :return: psi': list of wave function and wave-function derivative's first derivative values
    """
    dpsi = [0, 0]  # initialize the (psi', psi'') vector
    dpsi[0] = psi[1]  # set psi' = psi'
    if xi >= 0:
        dpsi[1] = ((xi+xi_naught) ** 2 - 2*epsilon) * psi[0]
    else:
        dpsi[1] = ((xi-xi_naught) ** 2 - 2*epsilon) * psi[0]
    return dpsi
