#!/usr/bin/env python3
# Gerald & Wheatley - Applied Numerical Analysis 7e
# Section 4.1 - Getting Chebyshev Polynomial
from sympy import Symbol, chebyshevt_poly, pprint 


def main():
    x = Symbol('x')

    for i in range(11):
        print(f"T_{i}(x):")
        pprint(chebyshevt_poly(i, x))
        print()


if __name__ == "__main__":
    main()
