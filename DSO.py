### DSO Eigenfunctions ###
from math import *
from numpy import *
import matplotlib.pyplot as mp

a0 = 1
a1 = 0
x0 = -0.5
xmin = 0
xmax = 4.0


def f(x):
    f = exp(-(x + x0) ** 2 / 2) * (a0 + a1 * (x + x0))  # analytical solution
    return f


x = linspace(xmin, xmax, 100)
mp.plot(x, f(x))


def f(x):
    f = exp(-(x - x0) ** 2 / 2) * (a0 + a1 * (x - x0))  # analytical solution
    return f


x = linspace(xmin, -xmax, 100)
mp.plot(x, f(x))
mp.show()
