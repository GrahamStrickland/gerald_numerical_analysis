#!/usr/bin/env python3
# Gerald & Wheatley - Applied Numerical Analysis 7e
# Section 4.2 - Pade Approximation
import numpy as np
from scipy.interpolate import pade
from sympy import Symbol, atan, pprint


def main():
    x = Symbol('x')

    maclaurin_series = atan(x).series(x, 0, 10)

    print("Maclaurin series approximation of arctan(x): ")
    pprint(maclaurin_series)

    coeffs = []
    for i in range(1, 10):
        coeffs.append(maclaurin_series.coeff(x**i).evalf())
    print(coeffs)

    p, q = pade(coeffs, 5)

    print(p)
    print(q)


if __name__ == "__main__":
    main()
