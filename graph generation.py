import matplotlib.pyplot as plt
import numpy as np

airfoil = 2412
data_file = open('Plunging_solution_results/result_file_NACA' + str(airfoil) + '.txt', 'r')
data_line = data_file.readlines()
aoa = float(data_line[5][:-1])
aoa = np.exp(1j * aoa)

file2_name = data_line[1][:-1] + '_data.txt'
airfoil_file = open(file2_name, 'r')
airfoil_line = airfoil_file.readlines()
airfoil_coordinate = airfoil_line[13][:-3].replace('i', 'j').split(' ')
airfoil_coordinate = np.array([complex(index) for index in airfoil_coordinate])
airfoil_coordinate = airfoil_coordinate * aoa

iteration = int(data_line[19][:-1])
current_iter_line = 25
time_step = float(data_line[15][:-1])

vortex_strength = data_line[27 + iteration * 2][1:-2].replace(' ', '').split(',')
vortex_strength = np.array([float(index) for index in vortex_strength])

pl_amplitude = float(data_line[7][:-1])
pl_frequency = float(data_line[9][:-1])
current_time = 0.0

for i in range(iteration):
    if i % 100 == 0:
        print('Iteration - ' + str(i))

    # plunging_dis = 1j * pl_amplitude * np.sin(2 * np.pi * pl_frequency * current_time)
    # plunging_vel = 2 * 1j * pl_amplitude * np.pi * pl_frequency * np.cos(2 * np.pi * pl_frequency * current_time)
    current_time += time_step
    data = data_line[current_iter_line][1:-2].replace(' ', '').split(',')
    data = np.array([complex(index) for index in data]) * aoa

    # plt.xlim(-5, 60)
    # plt.ylim(-20, 20)
    # plt.xlim(-0.1, 5)
    # plt.ylim(-5, 5)

    plt.axis('off')
    plt.grid(False)
    plt.plot(airfoil_coordinate.real, airfoil_coordinate.imag, color='b')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(data.real, data.imag, s=2, color='g')

    plt.savefig('Plunging_solution_results/' + str(i))
    plt.close()
    current_iter_line += 2
