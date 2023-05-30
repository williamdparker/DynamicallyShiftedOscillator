from wag_the_dog import wag_the_dog
from discrete_matrix import solve_discrete_eigenproblem
import numpy as np

# xi_0
shift = 0
# even or odd
expected_symmetry = 'odd'
# expected_symmetry = 'even'
# n
# eigenvalue_index = 1.0
# ε_expected
# expected_eigenvalue = eigenvalue_index + 0.5
# expected_eigenvalue = 0.309468
# expected_eigenvalue = 0.73422037
# expected_eigenvalue = 1.49997346
# expected_eigenvalue = 2.19740304
# expected_eigenvalue = 2.99866356
#
expected_eigenvalue = 1.4959
# expected_eigenvalue = 3.037

# δ
search_step = 1.e-4
search_starting_index, search_ending_index = -5, 5
# ε = ε_expected ± i δ
epsilon_values = expected_eigenvalue + np.arange(search_starting_index, search_ending_index+1, 1)*search_step

if expected_symmetry == 'even':
    # even psi0, psi0'
    initial_conditions = (1, 0)
else:
    # odd psi0, psi0'
    initial_conditions = (0, 1)

# [x_min, x_max]
#minimum_x, maximum_x = -5, 4
minimum_x, maximum_x = -10, 10

# N
# number_of_points = 2_000
# number_of_points = 5

# print(f'    N   ε_[0-4] (ħω)')
# for number_of_points in [100]:
for number_of_points in [1_000]:
# for number_of_points in [100, 200, 500, 1_000, 2_000, 5_000]:
# for number_of_points in [1_000, 2_000, 5_000]:
    positions = np.linspace(minimum_x, maximum_x, num=number_of_points)
#   wag_the_dog('dynamically shifted oscillator', epsilon_values, 1, positions, initial_conditions,
    solve_discrete_eigenproblem('dynamically shifted oscillator', 1, positions, reduced_action_quantum=1, mass=1,
                                starting_eigenvalue=0, number_of_eigenvalues=1, additional_parameters=[shift])

