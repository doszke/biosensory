import math
import numpy as np
from matplotlib import pyplot as plt


class Lambda:

    def __init__(self, km, k_cat, c_e0):
        self.km = km
        self.k_cat = k_cat
        self.c_e0 = c_e0

    def value(self):
        return self.k_cat*self.c_e0/self.km


class Equation:

    def __init__(self, n, delta_h, lmbda, s, r_eff, v):
        self.n = n
        self.delta_h = delta_h
        self.lmbda = lmbda
        self.s = s
        self.r_eff = r_eff
        self.v = v

    def calculate_c(self, delta_u, t):
        denominator = self.n * self.delta_h * self.lmbda.value() * self.s * self.r_eff * self.v * np.exp(-self.lmbda.value() * t)
        return delta_u / denominator

    def calculate_delta_u(self, c, t):
        denominator = self.n * self.delta_h * self.lmbda.value() * self.s * self.r_eff * self.v * np.exp(
            -self.lmbda.value() * t)
        return c * denominator


if __name__ == '__main__':
    lmbda = Lambda(km=1.1*10**-6, k_cat=8.0/150.0*10**-15, c_e0=1)
    eq = Equation(s=35*10**-6, r_eff=666.66, delta_h=-30*10**3, lmbda=lmbda, v=1000*10**-6, n=200)
    print(lmbda.value())
    time = np.zeros([100]) # 100 miejsc
    rt = 5 # time(100)
    dt = rt/100 # dt
    for x in range(100):
        time[x] = dt*x
    print(time)
    c = np.zeros([6])
    for x in range(6):
        c[x] = 49*10**-6 + x * 1*10**-6
    for x in range(6):
        data = np.zeros([100])
        for y in range(100):            # 6 serii po 100 pkt
            data[y] = eq.calculate_delta_u(c[x], time[y])
        plt.plot(time, data)
        print(data[99]-data[0])
    plt.show()
