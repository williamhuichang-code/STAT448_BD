# from Linear combinations, span, and basis vectors | Chapter 2, Essence of linear algebra
# https://www.youtube.com/watch?v=k7RM-ot2NWY&list=PLZHQObOWTQDPD3MizzM2xVFitgF8hE_ab&index=2

import numpy as np
import matplotlib.pyplot as plt

# Define two vectors in 2D
v = np.array([2, 1])
w = np.array([1, 3])

# Create scalar values for a
a_values = np.linspace(-5, 5, 100)

# Compute linear combinations: a*v + w
points = np.array([a * v + w for a in a_values])

# Plot
plt.figure(figsize=(6,6))

# Plot the resulting points
plt.plot(points[:,0], points[:,1], 'm-', label='a*v + w (varying a)')

# Plot vectors v and w
plt.quiver(0, 0, v[0], v[1], angles='xy', scale_units='xy', scale=1, color='red', label='v')
plt.quiver(0, 0, w[0], w[1], angles='xy', scale_units='xy', scale=1, color='blue', label='w')

plt.axhline(0)
plt.axvline(0)
plt.grid()
plt.legend()
plt.gca().set_aspect('equal', 'box')
plt.title("Tip of a*v + w traces a straight line")

plt.show()