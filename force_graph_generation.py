import matplotlib.pyplot as plt
import numpy as np

airfoil = 2412
# data_file = open('Plunging_solution_results/result_file_NACA' + str(airfoil) + '.txt', 'r')
data_file = open('Transient_solution_results/force_file_NACA' + str(airfoil) + '.txt', 'r')
data_line = data_file.readlines()
iteration = int(data_line[19])
time_step = float(data_line[15])

force_x = np.array([])
force_y = np.array([])
for i in range(iteration):
    a = data_line[i + 25][:-1].split()
    force_x = np.append(force_x, [a[0]])
    force_y = np.append(force_y, [a[1]])

    if i == 1:
        force_x[0] = a[0]
        force_y[0] = a[1]

time_list = time_step * np.arange(0, len(force_x))
# plt.plot(time_list, force_x)
plt.scatter(time_list, force_y)
plt.axis(enable=False)
plt.show()
