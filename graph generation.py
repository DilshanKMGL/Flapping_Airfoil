import matplotlib.pyplot as plt
import numpy as np

data_file = open('result_file.txt', 'r')
data_line = data_file.readlines()
aoa = float(data_line[5][:-1])
aoa = np.exp(-1j * aoa)

file2_name = data_line[1][:-1] + '_data.txt'
airfoil_file = open(file2_name, 'r')
airfoil_line = airfoil_file.readlines()
airfoil_coordinate = airfoil_line[11][:-3].replace('i', 'j').split(' ')
airfoil_coordinate = np.array([complex(index) for index in airfoil_coordinate])
airfoil_coordinate = airfoil_coordinate * aoa

iteration = int(data_line[19][:-1])
current_iter_line = 25

vortex_strength = data_line[27 + iteration * 2][1:-2].replace(' ', '').split(',')
vortex_strength = np.array([float(index) for index in vortex_strength])

for i in range(iteration):
    print('Iteration: ' + str(i+1))
    data = data_line[current_iter_line][1:-2].replace(' ', '').split(',')
    data = np.array([complex(index) for index in data]) * aoa

    plt.axis('off')
    plt.grid(False)
    plt.plot(airfoil_coordinate.real, airfoil_coordinate.imag, color='b')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(data.real, data.imag, s=2, color='g')

    plt.savefig(str(i))
    plt.close()
    current_iter_line += 2
