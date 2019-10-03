import numpy as np
import time
from matplotlib import pyplot as plt


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


def make_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
              time_step, current_time, iteration, distance, angle):
    heading = 'result_file.txt'
    file1 = open(heading, 'w')

    file1.write('airfoil\n' + str(airfoil) + '\n')
    file1.write('free_velocity\n' + str(free_velocity) + '\n')
    file1.write('free_aoa\n' + str(free_aoa) + '\n')
    file1.write('pl_amplitude\n' + str(pl_amplitude) + '\n')
    file1.write('pl_frequency\n' + str(pl_frequency) + '\n')
    file1.write('pi_amplitude\n' + str(pi_amplitude) + '\n')
    file1.write('pi_frequency\n' + str(pi_frequency) + '\n')
    file1.write('time_step\n' + str(time_step) + '\n')
    file1.write('current_time\n' + str(current_time) + '\n')
    file1.write('iteration\n' + str(iteration) + '\n')
    file1.write('distance\n' + str(distance) + '\n')
    file1.write('angle\n' + str(angle) + '\n')

    file1.close()


def update_file(te_vortex_z, iteration):
    heading = 'result_file.txt'
    file1 = open(heading, "a+")
    file1.write('te_vortex_z - iteration ' + str(iteration + 1) + '\n' + str(list(te_vortex_z)) + '\n')


def write_array(circulation, te_vortex_strength, iterate_time_step):
    heading = 'result_file.txt'
    file1 = open(heading, "a+")
    file1.write('circulation\n' + str(list(circulation)) + '\n')
    file1.write('te_vortex_strength\n' + str(list(te_vortex_strength)) + '\n')
    file1.write('iteration step time\n' + str(list(iterate_time_step)) + '\n')
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


start = time.time()
iterate_time = start
# ------ airfoil data
airfoil = 'NACA2412'
N, radius, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane = read_data(airfoil)
# ------ free stream velocity
re_num = 1.5e5
density = 1.225
viscosity = 1.789e-5
free_velocity = re_num * viscosity / density
free_aoa = 5.0  # in degrees

# ------ plunging parameters
pl_amplitude = 0
pl_frequency = 0
# ------ pitching parameters
pi_amplitude = 0
pi_frequency = 0
# ------ new vortex
distance = 0.01
angle = 0
angle = np.deg2rad(angle)
# ------ data store
circulation_list = np.array([])
te_vortex_strength = np.array([])
te_vortex_z = np.array([])
te_vortex_v = np.array([])
te_vortex_u = np.array([])

le_vortex_strength = np.array([])
le_vortex_z = np.array([])
le_vortex_v = np.array([])
le_vortex_u = np.array([])

iterate_time_step = np.array([])

time_step = 0.01  # " if the time step > 0.001, sudden variation of vortex position"
current_time = 0.00
iteration = 1

v_crit = 2.0  # LESP - critical velocity initialization
cal_les = True

# ----- write in a file
make_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
          time_step, current_time, iteration, distance, angle)

# ------ iteration code
for iterate in range(iteration):
    print('Iteration - ' + str(iterate + 1))

    # --- calculate velocity of the freestream - this will modify for plunging and flapping motion
    velocity = free_velocity
    aoa = np.deg2rad(free_aoa)

    # ------ calculate trailing edge position
    search_point = center_circle + radius
    trailing_edge_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, trailing_edge_z))
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)

    # ------ calculate leading edge position
    leading_edge_z = complex(min(z_plane))
    search_point = center_circle - radius
    leading_edge_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, leading_edge_z))
    leading_edge_u = complex((leading_edge_v - center_circle) / radius)

    # ------ calculate trailing edge vortex position
    new_te_vortex_position_z = trailing_edge_z + distance * pow(np.e, -1j * angle)
    new_te_vortex_position_z = complex(new_te_vortex_position_z)
    search_point = center_circle + radius
    new_te_vortex_position_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle,
                                              new_te_vortex_position_z))
    new_te_vortex_position_u = complex((new_te_vortex_position_v - center_circle) / radius)

    # ------ calculate trailing edge vortex strength and circulation
    # freestream component
    p1 = velocity * (np.exp(-1j * aoa) - np.exp(1j * aoa) / trailing_edge_u ** 2)
    # circulation component
    p2 = -1j / (2 * np.pi * trailing_edge_u)
    # newly sheded vortex
    p3 = -1j * (1 / (trailing_edge_u - new_te_vortex_position_u) +
                1 / (trailing_edge_u * (1 - trailing_edge_u * new_te_vortex_position_u))) / (2 * np.pi)
    # previously shed trailing edge vortices
    p4 = -1j * te_vortex_strength * (1 / (trailing_edge_u - te_vortex_u) +
                                     1 / (trailing_edge_u * (1 - trailing_edge_u * te_vortex_u))) / (2 * np.pi)
    p4 = sum(p4)  # convert array --> single_value
    # previously shed leading edge vortices
    p5 = -1j * le_vortex_strength * (1 / (trailing_edge_u - le_vortex_u) +
                                     1 / (trailing_edge_u * (1 - trailing_edge_u * le_vortex_u))) / (2 * np.pi)
    p5 = sum(p5)  # convert array --> single value
    # total vortex strength
    s = sum(te_vortex_strength) + sum(le_vortex_strength)
    strength_t = (p2 * s - p4 - p5 - p1) / (p3 - p2)
    circulation = -strength_t - s

    # ------ calculate leading edge velocity
    # - calculate complex_potential
    # freestream component
    p1 = velocity * (np.exp(-1j * aoa) - np.exp(1j * aoa) / leading_edge_u ** 2)
    # circulation component
    p2 = -1j * circulation / (2 * np.pi * leading_edge_u)
    # newly sheded vortex
    p3 = -1j * strength_t * (1 / (leading_edge_u - new_te_vortex_position_u) +
                             1 / (leading_edge_u * (1 - leading_edge_u * new_te_vortex_position_u))) / (2 * np.pi)
    # previously shed trailing edge vortices
    p4 = -1j * te_vortex_strength * (1 / (leading_edge_u - te_vortex_u) +
                                     1 / (leading_edge_u * (1 - leading_edge_u * te_vortex_u))) / (2 * np.pi)
    p4 = sum(p4)  # convert array --> single_value
    # previously shed leading edge vortices
    p5 = -1j * le_vortex_strength * (1 / (leading_edge_u - le_vortex_u) +
                                     1 / (leading_edge_u * (1 - leading_edge_u * le_vortex_u))) / (2 * np.pi)
    p5 = sum(p5)  # convert array --> single value

    # - calculate derivatives
    dwdu = p1 + p2 + p3 + p4 + p5
    power = np.arange(1, len(Gkn) + 1)
    dzdv = 1 - sum(power * Gkn * radius ** power / (leading_edge_v - center_circle) ** (power + 1))
    dudv = 1 / radius
    dudz = dudv / dzdv

    # Calculate leading edge velocity
    le_vel_conj = dwdu * dudz
    le_vel = np.conj(le_vel_conj)

    # ------ calculate leading edge vortex position
    new_le_vortex_position_z = leading_edge_z - distance * pow(np.e, -1j * angle)
    new_le_vortex_position_z = complex(new_le_vortex_position_z)
    search_point = center_circle - radius
    new_le_vortex_position_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle,
                                              new_le_vortex_position_z))
    new_le_vortex_position_u = complex((new_le_vortex_position_v - center_circle) / radius)

    print(abs(le_vel))
    if abs(le_vel) > v_crit and cal_les:
        # calculate trailing edge boundary condition
        # freestream component
        p1 = velocity * (np.exp(-1j * aoa) - np.exp(1j * aoa) / trailing_edge_u ** 2)
        # circulation component
        p2 = -1j * circulation / (2 * np.pi * trailing_edge_u)
        # newly shed vortex - trailing edge
        p3 = -1j * (1 / (trailing_edge_u - new_te_vortex_position_u) +
                    1 / (trailing_edge_u * (1 - trailing_edge_u * new_te_vortex_position_u))) / (2 * np.pi)
        # newly shed vortex - leading edge
        p4 = -1j * (1 / (trailing_edge_u - new_le_vortex_position_u) +
                    1 / (trailing_edge_u * (1 - trailing_edge_u * new_le_vortex_position_u))) / (2 * np.pi)
        # previously shed trailing edge vortices
        p5 = -1j * te_vortex_strength * (1 / (trailing_edge_u - te_vortex_u) +
                                         1 / (trailing_edge_u * (1 - trailing_edge_u * te_vortex_u))) / (2 * np.pi)
        p5 = sum(p5)  # convert array --> single_value
        # previously shed leading edge vortices
        p6 = -1j * le_vortex_strength * (1 / (trailing_edge_u - le_vortex_u) +
                                         1 / (trailing_edge_u * (1 - trailing_edge_u * le_vortex_u))) / (2 * np.pi)
        p6 = sum(p6)  # convert array --> single_value

        # calculate leading edge boundary condition
        # freestream component
        e1 = velocity * (np.exp(-1j * aoa) - np.exp(1j * aoa) / leading_edge_u ** 2)
        # circulation component
        e2 = -1j * circulation / (2 * np.pi * leading_edge_u)
        # newly shed vortex - trailing edge
        e3 = -1j * (1 / (leading_edge_u - new_te_vortex_position_u) +
                    1 / (leading_edge_u * (1 - leading_edge_u * new_te_vortex_position_u))) / (2 * np.pi)
        # newly shed vortex - leading edge
        e4 = -1j * (1 / (leading_edge_u - new_le_vortex_position_u) +
                    1 / (leading_edge_u * (1 - leading_edge_u * new_le_vortex_position_u))) / (2 * np.pi)
        # previously shed trailing edge vortices
        e5 = -1j * te_vortex_strength * (1 / (leading_edge_u - te_vortex_u) +
                                         1 / (leading_edge_u * (1 - leading_edge_u * te_vortex_u))) / (2 * np.pi)
        e5 = sum(e5)  # convert array --> single_value
        # previously shed leading edge vortices
        e6 = -1j * le_vortex_strength * (1 / (leading_edge_u - le_vortex_u) +
                                         1 / (leading_edge_u * (1 - leading_edge_u * le_vortex_u))) / (2 * np.pi)
        e6 = sum(e6)  # convert array --> single_value

        # equation set
        B1 = -(p5 + p6 + p1)
        B2 = np.conj(v_crit) / dudz - (e5 + e6 + e1)
        B3 = - sum(te_vortex_strength) - sum(le_vortex_strength)

        C1 = (e4 - e2) / (e3 - e2)
        C2 = (p4 - p2) / (p3 - p2)
        C3 = (B2 - e2 * B3) / (e3 - e2)
        C4 = (B1 - p2 * B3) / (p3 - p2)

        strength_l = (C3 - C4) / (C1 - C2)
        strength_t = C4 - strength_l * C2
        circulation = B3 - strength_l - strength_t

        le_vortex_z = np.append(le_vortex_z, [new_le_vortex_position_z])
        le_vortex_v = np.append(le_vortex_v, [new_le_vortex_position_v])
        le_vortex_u = np.append(le_vortex_u, [new_le_vortex_position_u])
        le_vortex_strength = np.append(le_vortex_strength, [strength_l])

    # ------ append all variables to corresponding arrays
    te_vortex_z = np.append(te_vortex_z, [new_te_vortex_position_z])
    te_vortex_v = np.append(te_vortex_v, [new_te_vortex_position_v])
    te_vortex_u = np.append(te_vortex_u, [new_te_vortex_position_u])
    te_vortex_strength = np.append(te_vortex_strength, [strength_t])
    circulation_list = np.append(circulation_list, [circulation])

    current_time += time_step

write_array(circulation_list, te_vortex_strength, iterate_time_step)
print('total time ', time.time() - start)

# plot trailing edge and leading edge vortices
plt.plot(z_plane.real, z_plane.imag)
plt.scatter(te_vortex_z.real, te_vortex_z.imag)
plt.scatter(le_vortex_z.real, le_vortex_z.imag)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

# plot selected leading and trailing edges of the airfoil and respective plane
'''
plt.plot(z_plane.real, z_plane.imag, color='b')
plt.scatter(trailing_edge_z.real, trailing_edge_z.imag, color='b')
plt.scatter(leading_edge_z.real, leading_edge_z.imag, color='b')
plt.scatter(new_le_vortex_position_z.real, new_le_vortex_position_z.imag, color='black')
plt.scatter(new_te_vortex_position_z.real, new_te_vortex_position_z.imag, color='black')

plt.plot(v_plane.real, v_plane.imag, color='g')
plt.scatter(trailing_edge_v.real, trailing_edge_v.imag, color='g')
plt.scatter(leading_edge_v.real, leading_edge_v.imag, color='g')
plt.scatter(new_le_vortex_position_v.real, new_le_vortex_position_v.imag, color='black')
plt.scatter(new_te_vortex_position_v.real, new_te_vortex_position_v.imag, color='black')

plt.plot(u_plane.real, u_plane.imag, color='r')
plt.scatter(trailing_edge_u.real, trailing_edge_u.imag, color='r')
plt.scatter(leading_edge_u.real, leading_edge_u.imag, color='r')
plt.scatter(new_le_vortex_position_u.real, new_le_vortex_position_u.imag, color='black')
plt.scatter(new_te_vortex_position_u.real, new_te_vortex_position_u.imag, color='black')

plt.gca().set_aspect('equal', adjustable='box')
plt.show()
'''
