import numpy as np

if __name__ == '__main__':
    # ξ0
    shift = -1

    # Position range
    # [ξ_min, ξ_max]
    minimum_xi, maximum_xi = -10, 10

    # Discretization
    number_of_points = 10_000
    delta_xi = 1 / number_of_points
    positions = np.linspace(minimum_xi, maximum_xi, num=number_of_points)

    # Unknown eigenvalues
    from discrete_matrix.discrete_matrix import solve_discrete_eigenproblem
    solve_discrete_eigenproblem('dynamically shifted oscillator', 1, positions, reduced_action_quantum=1, mass=1,
                                starting_eigenvalue=0, number_of_eigenvalues=5, additional_parameters=[shift],
                                write_file=False)

    # Guess eigenvalues
    ## Expected symmetry
    # symmetry = 'odd'
    # symmetry = 'even'
    # Set symmetry to pass to method
    # if expected_symmetry == 'even':
    #     # even psi0, psi0'
    #     initial_conditions = (1, 0)
    # else:
    #     # odd psi0, psi0'
    #     initial_conditions = (0, 1)


    ## RK4
    # from runge_kutta.wag_the_dog import wag_the_dog
    # wag_the_dog('dynamically shifted oscillator', epsilon_values, 1, positions, initial_conditions,
    #            additional_parameters=[shift])

    # EC
    # from euler_cromer.euler_cromer import euler_cromer



    ## Eigenvalue guess
    # expected_eigenvalue = eigenvalue_index + 0.5
    # expected_eigenvalue = 0.309468
    # expected_eigenvalue = 0.73422037
    # expected_eigenvalue = 1.4959
    # expected_eigenvalue = 1.496
    # expected_eigenvalue = 1.49997346
    # expected_eigenvalue = 2.19740304
    # expected_eigenvalue = 2.99866356
    # expected_eigenvalue = 3.037

    ## Range around guess
    # δ
    # search_step = 1.e-3
    # search_starting_index, search_ending_index = -5, 5
    # ε = ε_expected ± i δ
    # epsilon_values = expected_eigenvalue + np.arange(search_starting_index, search_ending_index+1, 1)*search_step





