import numpy as np
from matplotlib import pyplot as plt

radius = 0.27546
center_circle = 0.48333+0.01141j
num_div = 1000
phi = 2 * np.pi / num_div  # divisable number
angle = np.array([n*phi for n in range(num_div)])
circle_point = center_circle + radius * np.exp(1j * angle)
print(circle_point)
plt.plot(circle_point.real, circle_point.imag)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
