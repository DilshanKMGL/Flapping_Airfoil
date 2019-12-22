from matplotlib import pyplot as plt
import numpy as np

x1 = x2 = np.linspace(0,4*np.pi, 100)
y1 = np.sin(x1)
y2 = np.cos(x2)

plt.figure(1)
plt.plot(x1, y1, 'b-*')
plt.figure(2)
plt.plot(x2, y2, 'r-*')
plt.show()