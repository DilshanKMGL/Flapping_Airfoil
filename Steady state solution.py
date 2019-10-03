import numpy as np
import time


def make_file(airfoil, free_velocity, free_aoa, end_aoa, heading):
    file1 = open(heading, 'w')
    file1.write('airfoil\n' + str(airfoil) + '\n')
    file1.write('free_velocity\n' + str(free_velocity) + '\n')
    file1.write('free_aoa\n' + str(free_aoa) + '\n')
    file1.write('end_aoa\n' + str(end_aoa) + '\n')
    file1.write('angle\tcirculation\n')
    file1.close()


def read_data(heading):
    heading = 'Airfoil data/' + str(heading) + '_data.txt'
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


def write_array(angle, circulation, heading):
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
airfoil_list = ['0006', '0008', '0009', '0010', '0012', '0015', '0018', '0021', '0024', '1408', '1410', '1412', '2408',
                '2410', '2411', '2412', '2414', '2415', '2418', '2421', '2424', '4412', '4415', '4418', '4421', '4424',
                '6409', '6412']

re_num = 10e7
density = 1.225
viscosity = 1.789e-5
free_velocity = re_num * viscosity / density
free_aoa = -10
end_aoa = 19

for airfoil_number in airfoil_list:
    airfoil = 'NACA' + airfoil_number
    re_number_str = '10,'+str((len(str(re_num))-4))
    heading = 'Steady state solution results/' + airfoil + ' - ' + re_number_str + ' - steady state resuts.txt'
    N, radius, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane = read_data(airfoil)
    # ------ free stream velocity

    make_file(airfoil, free_velocity, free_aoa, end_aoa, heading)
    print(airfoil)

    for angle in np.arange(free_aoa, end_aoa, 0.1):
        # print('Angle of Attack ' + str(round(angle, 2)))
        angle_rad = np.deg2rad(angle)
        search_point = center_circle + radius
        trailing_edge_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, trailing_edge_z))
        trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)
        u1 = free_velocity * ((1 / pow(np.e, 1j * angle_rad)) - (pow(np.e, 1j * angle_rad) / (trailing_edge_u ** 2)))
        u2 = -1j / (2 * np.pi * trailing_edge_u)
        steady_circulation = complex(- u1 / u2)
        write_array(round(angle, 2), steady_circulation.real, heading)
