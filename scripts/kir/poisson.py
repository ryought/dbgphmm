#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
from scipy.special import factorial, comb
import matplotlib.pyplot as plt


def poisson(x, l):
    """
    average coverage is l.
    probability that there is k reads at a base.
    """
    a = np.power(l, x)
    b = np.exp(-l)
    c = factorial(x)
    return a * b / c


def error(p, k):
    return np.power(1-p, k)


def binom(c, k, p, x):
    """
    a k-mer is sequenced c times.
    the probability that there x k-mers sequenced without error.

    cCx (1-p)^{k*x}
    """
    e = error(p, k)
    return comb(c, x) * np.power(e, x) * np.power(1-e, c-x)
    # return comb(c, x) * np.power(1 - p, k * x)
    # return np.power(1 - p, k * x)


def kmer(c, k, p, t):
    """
    less than t times
    """
    return sum((binom(c, k, p, x) for x in range(0, t)))


def z(C, k, p):
    return sum((poisson(x, C) * kmer(x, k, p, 2) for x in range(1, 20)))


def main():
    C = 10
    ks = [k for k in range(10, 200)]
    p = 0.001
    G = 40000  # 100kb
    # zs = [z(C, k, p) * G for k in ks]
    plt.plot(ks, [error(0.001, k) for k in ks], label='0.1%')
    plt.plot(ks, [error(0.003, k) for k in ks], label='0.3%')
    plt.plot(ks, [error(0.01, k) for k in ks], label='1%')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
