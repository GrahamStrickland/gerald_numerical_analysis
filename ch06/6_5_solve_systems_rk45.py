#!/usr/bin/env python3
# Gerald & Wheatley - Applied Numerical Analysis 7e
# Section 6.5 - Solution to system of first-order equations
# dx/dt = xy + t, x(0) = 1, dy/dt = ty + x, y(0) = -1
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from scipy.integrate import solve_ivp 

matplotlib.rcParams['text.usetex'] = True


def system(t: float, y: np.ndarray) -> np.ndarray:
    return np.array([y[0]*y[1] + t, t*y[1] + y[0]])


def main() -> None:
    x = np.linspace(-1., 1., 25)
    y = np.linspace(-1., 1., 25)
    X, Y = np.meshgrid(x, y)

    dx, dy = system(0., np.array([X, Y]))

    fig, ax = plt.subplots()
    ax.quiver(X, Y, dx, dy)
    ax.set_title("Dynamics for system of differential equations\n" +
                r"$\begin{array}{lr}\frac{dx}{dt} = xy + t,\\ \frac{dy}{dt} = ty + x\end{array}$")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")

    y0 = np.array([1., -1.])

    solution = solve_ivp(
        fun=system,
        t_span=[0., .4],
        y0=y0,
        method='RK45',
        t_eval=np.arange(0., .5, .1)
        )
    print(solution)

    ax.plot(y0[0], y0[1], "ko")
    ax.plot(solution.y[0, :], solution.y[1, :], "ro-")

    plt.show()
    fig.savefig("6_5_solve_system_rk45.pdf")


if __name__ == "__main__":
    main()
