#!/usr/bin/env python3
# Gerald & Wheatley - Applied Numerical Analysis 7e
# Section 4.1 - Economizing a Series
from sympy import Poly, Symbol, chebyshevt_poly, exp, factorial, pprint


def main():
    x = Symbol('x')

    taylor_series = exp(x).series(x, 0, 7)

    print("Taylor series approximation of e^x: ")
    pprint(taylor_series)

    p = taylor_series.as_poly() - ((chebyshevt_poly(6, x) / factorial(6)) / 2**5)

    print("\nEconomization of degree 5: ")
    pprint(p.as_expr())


if __name__ == "__main__":
    main()
