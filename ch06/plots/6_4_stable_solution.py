#!/usr/bin/env python3
# Gerald & Wheatley - Applied Numerical Analysis 7e
# Section 6.4 - Solution to dy/dx = f(x,y) = -2y + 2, y(0) = -1
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True

# Determine sample of interval of definition I = (-inf, inf)
x = np.linspace(-1., 4.)

# Define function for analytical solution
def y(x):
    return 1. - 2.*np.exp(-2.*x)

# Define points from tables 6.11 and 6.12
eq_6_18 = np.array([np.linspace(.2, 4., 39), 
                    [-.34502, -.09946, .09477, 0.26264, .38971, .50675, .58701, .67195,
                      .71823, .78466, .80437, .86291, .85920, .91923, .89151, .96262,
                      .90646, 1.00004, .90645, 1.03746, .89146, 1.08087, .85911,
                      1.13723, .80422, 1.21554, .71801, 1.32834, .58667, 1.49367,
                      .38920, 1.73799, .09401, 2.10038, -.34613, 2.63884, -1.00168,
                      3.43951, -1.97749]])
eq_6_19 = np.array([np.linspace(0.0, 4., 41), 
                    [-1., -.6, -.28, -.024, .1808, .34464, .47571, .58057, .66446, 
                     .73156, .78525, .8282, .86256, .89005, .91204, .92963, .94371, 
                     .95496, .96397, .97118, .97694, .98155, .98524, .98819, .99056, 
                     .99244, .99396, .99516, .99613, .99691, .99752, .99802, .99842,
                     .99873, .99899, .99919, .99935, .99948, .99958, .99967, .99973]])

# Plot function on iterval
fig, ax = plt.subplots()
solutions = ax.plot(x, y(x), eq_6_18[0], eq_6_18[1], 'o', eq_6_19[0], eq_6_19[1], 'x')

# Add title, labels, and gridlines
plt.title("Gerald \& Wheatley - Section 6.4 - Stability Considerations for\n" +
           r"$\frac{dy}{dx} = f(x,y) = -2y + 2, y(0) = -1$")
plt.xlabel("$x$")
plt.ylabel("$y$")
ax.legend([r"$y(x) = 1 - 2e^{-2x}$", r"$y_{n+1} = y_{n-1} + 2hf(x_n, y_n)$",
           r"$y_{n+1} = y_{n+1} + y_n + hy_n' + O(h^2) \qquad$ (Simple Euler Method)"])

# Plot asymptotes
tlim = ax.get_xlim()
plt.hlines(0, tlim[0], tlim[1], colors='black', linestyles='solid')
Xlim = ax.get_ylim()
plt.vlines(0, Xlim[0], Xlim[1], colors='black', linestyles='solid')
ax.grid(visible=True)

# Show plot and save image
plt.show()
fig.savefig("6_4_stable_solution.pdf")
