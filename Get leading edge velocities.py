import numpy as np
from matplotlib import pyplot as plt
import time


def make_file(airfoil, free_velocity, free_aoa, end_aoa):
    heading = 'initial_circulation.txt'
    file1 = open(heading, 'w')

    file1.write('airfoil\n' + str(airfoil) + '\n')
    file1.write('free_velocity\n' + str(free_velocity) + '\n')
    file1.write('free_aoa\n' + str(free_aoa) + '\n')
    file1.write('end_aoa\n' + str(end_aoa) + '\n')
    file1.write('angle\tcirculation\n')
    file1.close()


def read_data(heading):
    heading = str(heading) + '_data.txt'
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


def write_array(angle, circulation):
    heading = 'initial_circulation.txt'
    file1 = open(heading, "a+")
    file1.write(str(angle) + '\t' + str(circulation) + '\n')
    file1.close()


def newton(x0, epsilon, max_iter, Gkn, radius, center_circle, equal_val):
    xn = x0
    for n in range(0, max_iter):
        power = np.arange(1, len(Gkn) + 1)
        fxn = xn + sum(Gkn * (radius / (xn - center_circle)) ** power) - equal_val
        if abs(fxn) < epsilon:
            # print('Found solution after', n, 'iterations.')
            return xn
        Dfxn = 1 - sum(power * Gkn * ((radius ** power) / (xn - center_circle) ** (power + 1)))
        if Dfxn == 0:
            # print('Zero derivative. No solution found.')
            return None
        xn = (xn - fxn / Dfxn)
    print('Exceeded maximum iterations. No solution found.')
    return None


# ------ airfoil data
airfoil = 'NACA2412'
N, radius, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane = read_data(airfoil)
# ------ free stream velocity
re_num = 1e5
density = 1.225
viscosity = 1.789e-5
free_velocity = 10  # re_num*viscosity/density
free_aoa = -10
end_aoa = 14
velocity_array = []
x_axis = []
make_file(airfoil, free_velocity, free_aoa, end_aoa)

for angle in np.arange(free_aoa, end_aoa, 0.2):
    x_axis.append(angle)
    print('Angle of Attack ' + str(round(angle, 2)))
    free_aoa = np.deg2rad(angle)
    search_point = center_circle + radius
    trailing_edge_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, trailing_edge_z))
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)
    u1 = free_velocity * ((1 / pow(np.e, 1j * free_aoa)) - (pow(np.e, 1j * free_aoa) / (trailing_edge_u ** 2)))
    u2 = -1j / (2 * np.pi * trailing_edge_u)
    steady_circulation = complex(- u1 / u2)

    leading_edge_z = min(z_plane)
    search_point = center_circle - radius
    leading_edge_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, leading_edge_z))
    leading_edge_u = complex((leading_edge_v - center_circle) / radius)
    leading_edge_vel = -1j * steady_circulation / (2 * np.pi * leading_edge_u) + \
                       free_velocity * ((1 / pow(np.e, 1j * free_aoa)) -
                                        (pow(np.e, 1j * free_aoa) / (leading_edge_u ** 2)))
    velocity_array.append(abs(leading_edge_vel))

    write_array(round(angle, 2), steady_circulation.real)

heading = 'initial_circulation.txt'
file1 = open(heading, "a+")
file1.write('leading edge circulation\n' + str(list(velocity_array)))
file1.close()

plt.plot(x_axis, velocity_array)
plt.grid()
plt.title('leading_edge_vel vs angle of attack')
plt.xlabel('aoa')
plt.ylabel('velocity')
plt.show()
