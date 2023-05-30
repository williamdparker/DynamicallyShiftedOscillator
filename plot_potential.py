import numpy as np
import matplotlib.pyplot as plt


def calculate_potential(position, x0=0, factor=1):
    return np.piecewise(position, [position < 0, position >= 0],
                        [lambda x: 0.5*factor*(x + x0)**2,
                         lambda x: 0.5*factor*(x - x0)**2])


figure, axes = plt.subplots()
shifts = np.arange(-1, 2, 1)
minimum_position, maximum_position = -3, 3
number_of_positions = 100
positions = np.linspace(minimum_position, maximum_position, num=number_of_positions)
for shift in shifts:
    potential_values = calculate_potential(positions, x0=shift)
    plt.plot(positions, potential_values)

axes.set_aspect('equal')
plt.ylim([0, 6])
plt.show()
