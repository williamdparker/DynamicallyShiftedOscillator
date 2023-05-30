import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def wag_the_dog(potential_name, adjustable_coefficient_values, length_scale, positions, initial_values,
                additional_parameters=[]):
    # List of currently supported potentials
    supported_potentials = ['dynamically shifted oscillator']

    # Set up tuple for the range of positions to use in solve_ivp
    position_range = (np.min(positions), np.max(positions))

    figure, axes = plt.subplots()
    if potential_name == 'dynamically shifted oscillator':
        from dynamically_shifted_oscillator import dynamically_shifted_oscillator_differential_equation, \
            draw_dynamically_shifted_oscillator_potential

        dynamical_shift = additional_parameters[0]
        for adjustable_coefficient in adjustable_coefficient_values:
            # Assign epsilon to an list
            arguments = [adjustable_coefficient, dynamical_shift]
            # Solve y''(x) = f[x, y(x), y'(x), K] given (x_min, x_max), (y(x_min), y'(x_min)), x, and K
            solution = solve_ivp(dynamically_shifted_oscillator_differential_equation, position_range, initial_values,
                                 t_eval=positions, args=arguments)
            # Assign the normalized solution to psi to the y-axis
            #   A = 1 / sqrt( int( y(x)^2 ) dx )
            normalization_factor = 1. / np.sqrt(np.trapz(solution.y[0]**2, x=positions))
            # normalization_factor = 1. / np.linalg.norm(solution.y[0])
            print(f'A = {normalization_factor}')
            #   psi(x) = A * y(x)
            psi_numerical = normalization_factor * solution.y[0]
            # Plot normalized solution
            plt.plot(positions, psi_numerical, label=r'$\epsilon = {:.5f}$'.format(adjustable_coefficient))

        plt.legend()
        potential_positions = np.linspace(-3*length_scale, 3*length_scale)
        # potential_positions = positions
        draw_dynamically_shifted_oscillator_potential(potential_positions, dynamical_shift, potential_width=length_scale)

    else:
        print('No potential named {} supported'.format(potential_name))
        print('Try one from this list: {}'.format([potential for potential in supported_potentials]))
        return

    # Formatting for all plots
    #    Draw psi = 0 line
    plt.axhline(color='black')

    if len(adjustable_coefficient_values) > 5:
        axes.legend(loc='upper right', bbox_to_anchor=(1.4, 1))
        figure.tight_layout()
    #   Label axes
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\psi_{n}(x)$')
    #   Show plot
    plt.show()

    return
