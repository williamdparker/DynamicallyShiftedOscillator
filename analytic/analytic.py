import numpy as np
from scipy.special import hermite


def psi(xi, xi0, n, parity='odd'):
    if parity == 'odd':
        argument = xi + xi0
    else:
        argument = xi - xi0
    return np.exp(-0.5*argument**2) * hermite(n)(argument)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    x_values = np.linspace(-10, 10, num=100)
    for index in [0, 1, 2]:
        # print(f'n = {index}, n%2 = {index % 2}')
        if (index % 2) == 0:
            plt.plot(x_values, psi(x_values, -1, index, parity='even'))
        else:
            plt.plot(x_values, psi(x_values, +1, index))
    plt.show()

