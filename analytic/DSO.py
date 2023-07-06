### DSO Eigenfunctions ###
from math import *
from numpy import *
import matplotlib.pyplot as mp

a0 = 1
a1 = 1
x0 = 1
xmin = -10
xmax = 10


def f(x):
    f = exp(-(x + x0) ** 2 / 2) * (a0 + a1 * (x + x0))  # analytical solution
    return f


x = linspace(xmin, xmax, 1_000)
mp.plot(x, f(x))


def f(x):
    f = exp(-(x - x0) ** 2 / 2) * (a0 + a1 * (x - x0))  # analytical solution
    return f


x = linspace(xmin, -xmax, 1_000)
mp.plot(x, f(x))
mp.show()
