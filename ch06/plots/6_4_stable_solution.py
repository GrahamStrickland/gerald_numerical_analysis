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

# Plot function on iterval
fig, ax = plt.subplots()
solutions = ax.plot(x, y(x))

# Add title, labels, and gridlines
plt.title("Gerald \& Wheatley - Section 6.4 - Stability Considerations")
plt.xlabel("$x$")
plt.ylabel("$y$")
ax.legend([r"$y(x) = 1 - 2e^{-2x}$"])

# Plot asymptotes
tlim = ax.get_xlim()
plt.hlines(0, tlim[0], tlim[1], colors='black', linestyles='solid')
Xlim = ax.get_ylim()
plt.vlines(0, Xlim[0], Xlim[1], colors='black', linestyles='solid')
ax.grid(visible=True)

# Show plot and save image
plt.show()
fig.savefig("6_4_stable_solution.pdf")
