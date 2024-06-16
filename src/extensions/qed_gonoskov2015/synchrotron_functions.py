# This code provides computation of constants and tests for synchrotron functions
from math import *
from pipic.extensions import qed_gonoskov2015
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

F_1 = pi * 2 ** (5 / 3.0) / (sqrt(3) * gamma(1 / 3.0))
F_2 = sqrt(pi / 2.0)
G_1 = F_1 / 2.0
G_2 = F_2

print("F_1 =", F_1)
print("F_2 =", F_2)
print("G_1 =", G_1)
print("G_2 =", G_2)

print("K_5_3(x -> 0) =", 0.5 * gamma(2 / 3) * 2 ** (2 / 3))

x = 0.001
v = qed_gonoskov2015.synch_func_2(x=x)
print(v / (x ** (1 / 3)))


def plot_synch1():
    x = np.linspace(0.000001, 7, 1000)
    y = np.linspace(0.000001, 7, 1000)
    for i in range(1000):
        y[i] = qed_gonoskov2015.synch_func_1(x=x[i])
        # y[i] = qed_aeg.synch_func_2(x=x[i])

    plt.plot(x, y)
    plt.savefig("synch1.png")


plot_synch1()
