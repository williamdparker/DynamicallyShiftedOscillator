import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler


def calculate_potential(position, x0=0, factor=1):
    return np.piecewise(position, [position < 0, position >= 0],
                        [lambda x: 0.5*factor*(x + x0)**2,
                         lambda x: 0.5*factor*(x - x0)**2])


if __name__ == '__main__':
    custom_cycler = (cycler(color=['tab:red', 'tab:orange', 'tab:blue', 'tab:purple']) +
                     cycler(linestyle=['-', '--', ':', '-.']))
    plt.rc('axes', prop_cycle=custom_cycler)
    figure, axes = plt.subplots()
    # shifts = np.arange(-1, 2, 1)
    # print(f'xis = {shifts}')
    # shifts = [-1]
    shifts = [-2, -1, 1, 2]




    minimum_position, maximum_position = -4, 4
    number_of_positions = 10_000
    positions = np.linspace(minimum_position, maximum_position, num=number_of_positions)
    for shift in shifts:
        # print(f'xi0 = {shift}')
        potential_values = calculate_potential(positions, x0=shift)
        plt.plot(positions, potential_values, label=f'{shift:2}')

    axes.set_aspect('equal')
    plt.ylim([0, 4])
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$V$ ($\hbar\omega$)')
    plt.legend(title=r'$\xi_0$')
    plt.tight_layout()
    # plt.savefig('DSOPotentialShiftPlusOne.png')
    plt.show()
