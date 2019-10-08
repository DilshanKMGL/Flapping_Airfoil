import numpy as np


def read_data(heading):
    heading = 'Airfoil_data/' + str(heading) + '_data.txt'
    file1 = open(heading, 'r')
    line = file1.readlines()

    N = int(line[1])
    r = float(line[3])
    center_circle = complex(line[5].replace('i', 'j'))
    trailing_edge_z = complex(line[7].replace('i', 'j'))
    Gkn = [complex(index.replace('i', 'j')) for index in line[9][:len(line[9]) - 3].split(' ')]
    z_plane = np.asarray([complex(index.replace('i', 'j')) for index in line[11][:len(line[11]) - 3].split(' ')])
    v_plane = np.asarray([complex(index.replace('i', 'j')) for index in line[13][:len(line[13]) - 3].split(' ')])
    u_plane = np.subtract(v_plane, center_circle) / r

    return N, r, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane


N, radius, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane = read_data('NACA2412')

equal_val = trailing_edge_z
epsilon = 1e-8
xn = center_circle + radius
power = np.arange(1, len(Gkn) + 1)

while True:
    fxn = xn + sum(Gkn * (radius / (xn - center_circle)) ** power) - equal_val
    if abs(fxn) < epsilon:
        break
    dfxn = 1 - sum(Gkn * power * (radius ** power / (xn - center_circle)) ** (power + 1))
    xn -= fxn / dfxn
print(xn)
# (0.7589761712369115+0.002230175020331722j)