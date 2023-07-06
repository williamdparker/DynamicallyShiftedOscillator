from matplotlib import pyplot as plt


def euler_cromer(maximum_x, delta_x, x_naught, epsilon, parity,
                 convergence_threshold=1e-3):
    if parity == 'even':
        psi = 1
        dpsi = 0
    else:
        psi = 0
        dpsi = 1
    difference_from_maximum = 1
    x = 0
    x_s = [x]
    psis = [psi]
    while difference_from_maximum > convergence_threshold:
        d2psi = ((x + x_naught)*(x + x_naught) - 2*epsilon) * psi
        dpsi += d2psi * delta_x
        psi += dpsi * delta_x
        x += delta_x
        x_s.append(x)
        psis.append(psi)
        difference_from_maximum = maximum_x - x
        # print(difference_from_maximum)

    current_curve = plt.plot(x_s, psis, label=f'{epsilon:.5f}')
    if parity == 'even':
        plt.plot(-1*np.array(x_s), psis, color=current_curve[0].get_color())
    else:
        plt.plot(-1*np.array(x_s), -1*np.array(psis), color=current_curve[0].get_color())

    return


if __name__ == '__main__':
    import numpy as np
    import matplotlib
    matplotlib.use('macosx')
    test_epsilons = np.linspace(1.4997, 1.4998, 2)
    for test_epsilon in test_epsilons:
        euler_cromer(10, 0.001, -1.0, test_epsilon, 'even', convergence_threshold=1e-4)

    plt.xlim([-7, 7])
    plt.ylim([-1.5, 1.5])
    plt.legend(title=r'$\epsilon$')

    plt.show()



