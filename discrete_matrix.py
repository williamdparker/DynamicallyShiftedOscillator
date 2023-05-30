import numpy as np
import matplotlib.pyplot as plt


def solve_discrete_eigenproblem(potential_name, length_scale, positions, reduced_action_quantum=1, mass=1,
                                starting_eigenvalue=0, number_of_eigenvalues=3, additional_parameters=[]):

    # assuming equally spaced positions as grid
    number_of_grid_points = len(positions)
    grid_spacing = positions[1] - positions[0]
    kinetic_energy_scale = reduced_action_quantum**2 / (2 * mass * grid_spacing**2)

    upper_ones = np.diagflat(np.ones(number_of_grid_points-1), 1)
    lower_ones = np.diagflat(np.ones(number_of_grid_points-1), -1)

    if potential_name == 'dynamically shifted oscillator':
        force_constant = 1/length_scale
        dynamical_shift = additional_parameters[0]
        potential_grid = np.piecewise(positions, [positions < 0, positions >= 0],
                                      [lambda x: 0.5 * force_constant * (x + dynamical_shift) ** 2,
                                       lambda x: 0.5 * force_constant * (x - dynamical_shift) ** 2])

    diagonal_values = np.diagflat(2 + potential_grid/kinetic_energy_scale)
    hamiltonian_matrix = kinetic_energy_scale*(-upper_ones + diagonal_values - lower_ones)

    eigenvalues, eigenvectors = np.linalg.eig(hamiltonian_matrix)
    sorted_eigenvalue_indices = np.argsort(eigenvalues)
    selected_eigenvalue_indices = sorted_eigenvalue_indices[starting_eigenvalue:
                                                            starting_eigenvalue+number_of_eigenvalues]
    selected_eigenvalue_eigenvectors = eigenvectors[:, selected_eigenvalue_indices]
    selected_eigenvalue_eigenvectors = selected_eigenvalue_eigenvectors.transpose()

    for index, eigenvector in enumerate(selected_eigenvalue_eigenvectors):
        if np.sum(eigenvector) < 0:
            eigenvector = -1 * eigenvector
        eigenvector = eigenvector / np.linalg.norm(eigenvector)
        plt.plot(positions, eigenvector, label=f'{index+starting_eigenvalue}')

    axis_limits = [plt.xlim(), plt.ylim()]
    axis_ranges = np.diff(axis_limits, axis=1)
    text_position = [axis_limits[0][0] + 0.05*axis_ranges[0][0], axis_limits[1][0] + 0.9*axis_ranges[1][0]]
    plt.text(text_position[0], text_position[1], rf'$\xi_0 = {dynamical_shift}$')
    plt.legend(title=r'$n$')
    plt.show()
    print(f'{number_of_grid_points:6}\t{np.sort(eigenvalues)[:5]}')

    return eigenvalues, eigenvectors
