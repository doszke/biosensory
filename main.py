import numpy as np
from matplotlib import pyplot as plt


class Equation:

    def __init__(self, n, delta_h, lmbda, s, r_eff, v):
        self.n = n
        self.delta_h = delta_h
        self.lmbda = lmbda
        self.s = s
        self.r_eff = r_eff
        self.v = v

    def calculate_c(self, delta_u, t):
        denominator = self.n * self.delta_h * self.lmbda * self.s * self.r_eff * self.v * np.exp(-self.lmbda * t)
        return delta_u / denominator

    def calculate_delta_u(self, c, t):
        denominator = self.n * self.delta_h * self.lmbda * self.s * self.r_eff * self.v * np.exp(
            -self.lmbda * t)
        return c * denominator


def delta_u_from_lambda(eq):
    before_lmbda = eq.lmbda
    c = 50 * 10 ** -6  # 50 nmol/ ml
    base_lambda = 0.05
    t = 10
    lmbda = np.zeros([11])
    for x in range(11):
        lmbda[x] = 0.5 * base_lambda + 0.1 * base_lambda * x
    deltau = np.zeros([11])
    for x in range(11):
        eq.lmbda = lmbda[x]
        deltau[x] = eq.calculate_delta_u(c, t)
    plt.plot(lmbda, deltau)
    plt.xlabel("lambda [1/s]")
    plt.ylabel("ΔU [V]")
    eq.lmbda = before_lmbda
    plt.show()


def delta_u_from_c(eq):
    dt = 0.5
    t = np.zeros(11)
    for x in range(0, 11):
        t[x] = (x) * dt
    c = np.zeros([50])
    for y in range(0, 11):
        for x in range(50):
            c[x] = 1 * 10 ** -6 * x
        data = eq.calculate_delta_u(c, t[y])
        plt.plot(c, data)
    plt.xlabel("c [M]")
    plt.ylabel("ΔU [V]")
    plt.show()


def delta_u_from_r_eff(eq):
    before_reff = eq.r_eff
    c = 50 * 10 ** -6  # 50 nmol/ ml
    t = 10
    base_reff = 666.66
    r_eff = np.zeros([11])
    for x in range(11):
        r_eff[x] = 0.5 * base_reff + 0.1 * base_reff * x
    deltau = np.zeros([11])
    for x in range(11):
        eq.r_eff = r_eff[x]
        deltau[x] = eq.calculate_delta_u(c, t)
    plt.plot(r_eff, deltau)
    plt.xlabel("effective thermal resistance [K/W]")
    plt.ylabel("ΔU [V]")
    eq.r_eff = before_reff
    plt.show()


def delta_u_from_cm(eq):
    c = 50 * 10 ** -6  # 50 nmol/ ml
    t = 10
    c_vec = np.zeros([11])
    for x in range(11):
        c_vec[x] = 0.5 * c + 0.1 * c * x
    deltau = np.zeros([11])
    for x in range(11):
        deltau[x] = eq.calculate_delta_u(c_vec[x], t)
    plt.plot(c_vec, deltau)
    plt.xlabel("concentration [M]")
    plt.ylabel("ΔU [V]")
    plt.show()


if __name__ == '__main__':
    eq = Equation(s=35*10**-6, r_eff=666.66, delta_h=-30*10**3, lmbda=0.05, v=1000*10**-6, n=200)
    delta_u_from_cm(eq)
    delta_u_from_lambda(eq)
    delta_u_from_r_eff(eq)
