#!/usr/bin/env python3
# Gerald & Wheatley - Applied Numerical Analysis 7e
# Section 4.1 - Economizing a Series
from sympy import Symbol, chebyshevt_poly, exp, factorial, pprint


def main():
    x = Symbol('x')

    maclaurin_series = exp(x).series(x, 0, 7)

    print("Maclaurin series approximation of e^x: ")
    pprint(maclaurin_series)

    p = maclaurin_series.as_poly() - ((chebyshevt_poly(6, x) / factorial(6)) / 2**5)

    print("\nEconomization of degree 5: ")
    pprint(p.as_expr())


if __name__ == "__main__":
    main()
