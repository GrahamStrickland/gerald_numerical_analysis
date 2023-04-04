#!/usr/bin/env python3
# Gerald & Wheatley - Applied Numerical Analysis 7e
# Section 4.1 - Getting Chebyshev Polynomial
import numpy as np
from sympy import chebyshevt_poly, Symbol


def main():
    x = Symbol('x')

    for i in range(11):
        print(f"T_{i}(x) = {chebyshevt_poly(i, x)}")


if __name__ == "__main__":
    main()
