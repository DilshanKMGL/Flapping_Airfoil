import graph
from matplotlib import pyplot as plt
import numpy as np


y = []
z = []

iteration = 1000
x = np.linspace(0, 10*np.pi, iteration)
nrow = 2
ncol = 2
fig, axs = plt.subplots(ncols=ncol, nrows=nrow)
for i in range(iteration):
    y.append(np.sin(x[i]))
    z.append(np.cos(x[i]))
    islast = False
    if i == iteration - 1:
        islast = True
    temp_x = x[:i+1]
    graph.main(fig, axs, nrow, ncol, temp_x, y, temp_x, z, 0.005, islast)
    # graph.main(temp_x, y, temp_x, z, 0.005, islast)
