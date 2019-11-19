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

te_vortex_strength = data_line[27 + iteration * 2][1:-2].replace(' ', '').split(',')
te_vortex_strength = np.array([float(index) for index in te_vortex_strength])

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

    te_strength_part = te_vortex_strength[:len(data)]
    # print(te_strength_part)
    cw_vortex = np.array([])
    ccw_vortex = np.array([])

    for j in range(len(te_strength_part)):
        if te_strength_part[j] > 0:
            ccw_vortex = np.append(ccw_vortex, [data[j]])
        else:
            cw_vortex = np.append(cw_vortex, [data[j]])

    plt.xlim(-1, 20)
    plt.ylim(-8, 6)

    plt.axis('off')
    plt.grid(False)
    plt.plot(airfoil_coordinate.real, airfoil_coordinate.imag, color='g')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(ccw_vortex.real, ccw_vortex.imag, s=2, color='b')
    plt.scatter(cw_vortex.real, cw_vortex.imag, s=2, color='r')

    plt.savefig('Plunging_solution_results/' + str(i))
    plt.close()
    current_iter_line += 2
