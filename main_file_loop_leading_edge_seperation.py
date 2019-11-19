import numpy as np
import time


def read_data(heading):
    heading = 'Airfoil_data/' + str(heading) + '_data.txt'
    file1 = open(heading, 'r')
    line = file1.readlines()

    N = int(line[1])
    r = float(line[3])
    center_circle = complex(line[5].replace('i', 'j'))
    trailing_edge_z = complex(line[7].replace('i', 'j'))
    trailing_edge_v = complex(line[9].replace('i', 'j'))
    Gkn = np.asarray([complex(index.replace('i', 'j')) for index in line[11][:len(line[11]) - 3].split(' ')])
    z_plane = np.asarray([complex(index.replace('i', 'j')) for index in line[13][:len(line[13]) - 3].split(' ')])
    v_plane = np.asarray([complex(index.replace('i', 'j')) for index in line[15][:len(line[15]) - 3].split(' ')])
    u_plane = np.subtract(v_plane, center_circle) / r

    return N, r, center_circle, trailing_edge_z, trailing_edge_v, Gkn, z_plane, v_plane, u_plane


def make_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
              time_step, current_time, iteration, distance, angle, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
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


def update_file(te_vortex_z, le_vortex_z, iteration, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write('te_vortex_z - iteration ' + str(iteration + 1) + '\n' + str(list(te_vortex_z)) + '\n')
    file1.write('le_vortex_z - iteration ' + str(iteration + 1) + '\n' + str(list(le_vortex_z)) + '\n')


def write_array(circulation, te_vortex_strength, iterate_time_step, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write('circulation\n' + str(list(circulation)) + '\n')
    file1.write('te_vortex_strength\n' + str(list(te_vortex_strength)) + '\n')
    file1.write('iteration step time\n' + str(list(iterate_time_step)) + '\n')
    file1.close()


# mapping function should be validated
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
        xn -= fxn / Dfxn
    print('Exceeded maximum iterations. No solution found.')
    return None


# remove diagonal element from array and traspose the result
def diag_remove(p):
    """
    make square matrix from the given array and remove diagonal elements
    :param p: array
    :return: nd array
    """
    strided = np.lib.stride_tricks.as_strided

    p = np.tile(p, (len(p), 1))
    m = p.shape[0]
    s0, s1 = p.strides
    p = strided(p.ravel()[1:], shape=(m - 1, m), strides=(s0 + s1, s1)).reshape(m, -1)
    return p


start = time.time()
iterate_time = start
# ------ airfoil data
# 2410 2418
airfoil = 'NACA2412'
N, radius, center_circle, trailing_edge_z, trailing_edge_v, Gkn, z_plane, v_plane, u_plane = read_data(airfoil)
# ------ free stream velocity
re_num = 1e6
density = 1.225
viscosity = 1.789e-5
free_velocity = re_num * viscosity / density
free_aoa = 0.0
free_aoa = np.deg2rad(free_aoa)
# ------ plunging parameters
pl_amplitude = 0.5
pl_frequency = 5
# ------ pitching parameters
pi_amplitude = 0
pi_frequency = 0
# ------ new vortex
distance = 0.001
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

# activate several functions
plunging_on = True
leading_edge_emmit = True
leading_edge_velocity_crit = 0.06

iterate_time_step = np.array([])
# ------ time step
time_step = 0.001
current_time = 0.00
iteration = 1000

heading_file = 'Plunging_solution_results/' + 'result_file_' + airfoil + '.txt'
# ----- write in a file

make_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
          time_step, current_time, iteration, distance, angle, heading_file)
print(airfoil)

velocity_file = open('velocity.txt', 'w')

# ------ iteration code
for iterate in range(iteration):
    if iterate % 100 == 0:
        print('Iteration - ' + str(iterate))

    # ------ calculate trailing edge position and leading edge position
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)

    leading_edge_z = min(z_plane)
    leading_edge_v = newton(min(v_plane), 1e-8, 250, Gkn, radius, center_circle, leading_edge_z)
    leading_edge_u = complex((leading_edge_v - center_circle) / radius)

    if plunging_on:
        # ------ calculate plunging parameters
        plunging_dis = 1j * pl_amplitude * np.sin(2 * np.pi * pl_frequency * current_time)
        plunging_vel = 2 * 1j * pl_amplitude * np.pi * pl_frequency * np.cos(2 * np.pi * pl_frequency * current_time)

        # --- calculate velocity
        velocity = np.abs(free_velocity + plunging_vel)
        aoa = free_aoa + np.angle(free_velocity + plunging_vel)
    else:
        velocity = free_velocity
        aoa = free_aoa

    # ------ calculate new vortex position - trailing edge
    new_vortex_position_te_z = trailing_edge_z + distance * pow(np.e, -1j * angle)
    new_vortex_position_te_z = complex(new_vortex_position_te_z)
    search_point = center_circle + radius
    new_vortex_position_te_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle,
                                              new_vortex_position_te_z))
    new_vortex_position_te_u = complex((new_vortex_position_te_v - center_circle) / radius)

    # ------ calculate new vortex position - leading edge
    new_vortex_position_le_z = leading_edge_z - distance * pow(np.e, -1j * angle)
    new_vortex_position_le_z = complex(new_vortex_position_le_z)
    search_point = center_circle - radius
    new_vortex_position_le_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle,
                                              new_vortex_position_le_z))
    new_vortex_position_le_u = complex((new_vortex_position_le_v - center_circle) / radius)

    # ------------------------------------- #
    #    calculate leading edge velocity    #
    # ------------------------------------- #
    # freestream calculation
    d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / leading_edge_u ** 2)

    # circulation
    if len(circulation_list) > 0:
        d2 = -1j * circulation_list[-1] / (2 * np.pi * leading_edge_u)
    else:
        d2 = 0

    # previously shed vortices - leading edge
    if len(le_vortex_strength) > 0:
        p1 = 1.0 / (leading_edge_u - le_vortex_u)
        p2 = 1.0 / (leading_edge_u * (1.0 - leading_edge_u * np.conjugate(le_vortex_u)))
        d31 = -1j * le_vortex_strength / (2 * np.pi) * (p1 + p2)
        d31 = sum(d31)
    else:
        d31 = 0

    # previously shed vortices - trailing edge
    if len(le_vortex_strength) > 0:
        p1 = 1.0 / (leading_edge_u - te_vortex_u)
        p2 = 1.0 / (leading_edge_u * (1.0 - leading_edge_u * np.conjugate(te_vortex_u)))
        d32 = -1j * te_vortex_strength / (2 * np.pi) * (p1 + p2)
        d32 = sum(d32)
    else:
        d32 = 0
    dwdu = d1 + d2 + d31 + d32

    # calculate derivatives
    power = np.arange(1, len(Gkn) + 1)
    dzdu = radius - sum(Gkn * power / leading_edge_u ** (power + 1))
    leading_edge_velocity_z = dwdu / dzdu
    velocity_file.write(str(leading_edge_velocity_z) + str('\n'))

    # ------ create function to calculate circulation
    if leading_edge_velocity_z > leading_edge_velocity_crit and leading_edge_emmit:
        # ----- leading edge calculations
        # velocity calculation - leading edge point
        d1L = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / leading_edge_u ** 2)
        # circulation - leading edge point
        d2L = -1j / (2 * np.pi * leading_edge_u)
        # newly sheded vortex - leading edge - leading edge point
        p1 = 1.0 / (leading_edge_u - new_vortex_position_le_u)
        p2 = 1.0 / (leading_edge_u * (1 - leading_edge_u * np.conj(new_vortex_position_le_u)))
        d31L = -1j / (2 * np.pi) * (p1 + p2)
        # newly sheded vortex - trailing edge - leading edge point
        p1 = 1.0 / (leading_edge_u - new_vortex_position_te_u)
        p2 = 1.0 / (leading_edge_u * (1 - leading_edge_u * np.conj(new_vortex_position_te_u)))
        d32L = -1j / (2 * np.pi) * (p1 + p2)
        # previously shed vortices - leading edge
        p1 = 1.0 / (leading_edge_u - le_vortex_u)
        p2 = 1.0 / (leading_edge_u * (1.0 - leading_edge_u * np.conjugate(le_vortex_u)))
        d41L = -1j * le_vortex_strength / (2 * np.pi) * (p1 + p2)
        d41L = sum(d41L)
        # previously shed vortices - trailing edge
        p1 = 1.0 / (leading_edge_u - te_vortex_u)
        p2 = 1.0 / (leading_edge_u * (1.0 - leading_edge_u * np.conjugate(te_vortex_u)))
        d42L = -1j * te_vortex_strength / (2 * np.pi) * (p1 + p2)
        d42L = sum(d42L)

        # ----- trailing edge calculation
        # velocity calculation - trailing edge point
        d1T = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / trailing_edge_u ** 2)
        # circulation - trailing edge point
        d2T = -1j / (2 * np.pi * trailing_edge_u)
        # newly sheded vortex - leading edge - trailing edge point
        p1 = 1.0 / (trailing_edge_u - new_vortex_position_le_u)
        p2 = 1.0 / (trailing_edge_u * (1 - trailing_edge_u * np.conj(new_vortex_position_le_u)))
        d31T = -1j / (2 * np.pi) * (p1 + p2)
        # newly sheded vortex - trailing edge - trailing edge point
        p1 = 1.0 / (trailing_edge_u - new_vortex_position_te_u)
        p2 = 1.0 / (trailing_edge_u * (1 - trailing_edge_u * np.conj(new_vortex_position_te_u)))
        d32T = -1j / (2 * np.pi) * (p1 + p2)
        # previously shed vortices - leading edge
        p1 = 1.0 / (trailing_edge_u - le_vortex_u)
        p2 = 1.0 / (trailing_edge_u * (1.0 - trailing_edge_u * np.conjugate(le_vortex_u)))
        d41T = -1j * le_vortex_strength / (2 * np.pi) * (p1 + p2)
        d41T = sum(d41T)
        # previously shed vortices - trailing edge
        p1 = 1.0 / (trailing_edge_u - te_vortex_u)
        p2 = 1.0 / (trailing_edge_u * (1.0 - trailing_edge_u * np.conjugate(te_vortex_u)))
        d42T = -1j * te_vortex_strength / (2 * np.pi) * (p1 + p2)
        d42T = sum(d42T)

        E1 = -(d1T + d41T + d42T)
        E2 = 1 / dzdu
        E3 = leading_edge_velocity_crit / E2 - (d1L + d41L + d42L)
        E4 = -(sum(te_vortex_strength) + sum(le_vortex_strength))
        E5 = d31T / d2T - 1
        E6 = d32T / d2T - 1
        E7 = E1 / d2T - E4
        E8 = d31L / d2L - 1
        E9 = d32L / d2L - 1
        E10 = E3 / d2L - E4
        E11 = E6 / E5
        E12 = E7 / E5
        E13 = E9 / E8
        E14 = E10 / E8

        # - vortices calculation
        te_vortex_strength_new = (E14 - E12) / (E13 - E11)
        le_vortex_strength_new = E12 - E11 * te_vortex_strength_new
        circulation = E4 - te_vortex_strength_new - le_vortex_strength_new
        # - circulation append
        circulation_list = np.append(circulation_list, [circulation])
        # - trailing edge append
        te_vortex_strength = np.append(te_vortex_strength, [te_vortex_strength_new])
        te_vortex_z = np.append(te_vortex_z, [new_vortex_position_te_z])
        te_vortex_v = np.append(te_vortex_v, [new_vortex_position_te_v])
        te_vortex_u = np.append(te_vortex_u, [new_vortex_position_te_u])
        # - leading edge append
        le_vortex_strength = np.append(le_vortex_strength, [le_vortex_strength_new])
        le_vortex_z = np.append(le_vortex_z, [new_vortex_position_le_z])
        le_vortex_v = np.append(le_vortex_v, [new_vortex_position_le_v])
        le_vortex_u = np.append(le_vortex_u, [new_vortex_position_le_u])

    else:
        s = sum(te_vortex_strength)
        # velocity calculation
        d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / trailing_edge_u ** 2)
        # circulation
        d2 = -1j / (2 * np.pi * trailing_edge_u)
        # newly sheded vortex
        p1 = 1.0 / (trailing_edge_u - new_vortex_position_te_u)
        p2 = 1.0 / (trailing_edge_u * (1 - trailing_edge_u * np.conj(new_vortex_position_te_u)))
        d3 = -1j / (2 * np.pi) * (p1 + p2)

        # previously shed vortices
        p1 = 1.0 / (trailing_edge_u - te_vortex_u)
        p2 = 1.0 / (trailing_edge_u * (1.0 - trailing_edge_u * np.conjugate(te_vortex_u)))
        d4 = -1j * te_vortex_strength / (2 * np.pi) * (p1 + p2)
        d4 = sum(d4)

        circulation = complex((s * d3 - d1 - d4) / (d2 - d3)).real

        circulation_list = np.append(circulation_list, [circulation])
        te_vortex_strength = np.append(te_vortex_strength, [- s - circulation])
        te_vortex_z = np.append(te_vortex_z, [new_vortex_position_te_z])
        te_vortex_v = np.append(te_vortex_v, [new_vortex_position_te_v])
        te_vortex_u = np.append(te_vortex_u, [new_vortex_position_te_u])

    # ------------------------------------- #
    #        update time and file           #
    # ------------------------------------- #
    iterate_time_step = np.append(iterate_time_step, [time.time() - iterate_time])
    update_file(te_vortex_z, le_vortex_z, iterate, heading_file)

    # - calculate derivative for trailing edge vortices
    power = np.arange(1, len(Gkn) + 1)
    Gkn_coeff_te = np.tile(Gkn, (len(te_vortex_u), 1)).transpose()
    power_coeff_te = np.tile(power, (len(te_vortex_u), 1)).transpose()
    te_coeff = np.tile(te_vortex_u, (len(Gkn), 1))

    dudv_te = 1 / radius
    dzdv_te = 1.0 - sum(Gkn_coeff_te * power_coeff_te / radius / te_coeff ** (power_coeff_te + 1))
    dvdz_te = 1 / dzdv_te
    dudz_te = dudv_te * dvdz_te

    d2zdv2_te = sum(Gkn_coeff_te * power_coeff_te * (power_coeff_te + 1) / radius ** 2 /
                    te_coeff ** (power_coeff_te + 2))
    d2udz2_te = - dudv_te * d2zdv2_te / dzdv_te ** 3

    # - calculate derivative for leading edge vortices
    if len(le_vortex_strength) > 0:
        power = np.arange(1, len(Gkn) + 1)
        Gkn_coeff_le = np.tile(Gkn, (len(le_vortex_u), 1)).transpose()
        power_coeff_le = np.tile(power, (len(le_vortex_u), 1)).transpose()
        le_coeff = np.tile(le_vortex_u, (len(Gkn), 1))

        dudv_le = 1 / radius
        dzdv_le = 1.0 - sum(Gkn_coeff_le * power_coeff_le / radius / le_coeff ** (power_coeff_le + 1))
        dvdz_le = 1 / dzdv_le
        dudz_le = dudv_le * dvdz_le

        d2zdv2_le = sum(Gkn_coeff_le * power_coeff_le * (power_coeff_le + 1) / radius ** 2 /
                        le_coeff ** (power_coeff_le + 2))
        d2udz2_le = - dudv_le * d2zdv2_le / dzdv_le ** 3
    else:
        dudz_le = 0
        d2udz2_le = 0

    # ------------------------------------- #
    #           move vortices               #
    # ------------------------------------- #

    # - move vortices - trailing edge emmits
    # velocity
    d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / te_vortex_u ** 2)
    # circulation
    d2 = -1j * circulation / (2 * np.pi * te_vortex_u)
    # leading edge emits vortices
    if len(le_vortex_z) > 0:
        le_ss = np.tile(le_vortex_strength, (len(le_vortex_strength), 1)).transpose()
        le_uu = np.conjugate(np.tile(le_vortex_u, (len(le_vortex_u), 1)).transpose())
        le_u = np.tile(le_vortex_u, (len(le_vortex_u), 1))
        d4L = sum(-1j * le_ss / (2 * np.pi) * (1 / (le_u * (le_u - le_uu)))) + \
              sum(-1j * le_ss / (2 * np.pi) * (1 / (le_u - le_uu)))
    else:
        d4L = 0

    # trailing edge emits vortices
    te_ss = diag_remove(te_vortex_strength).transpose()  # strength of vortices, remove diagonal and transpose
    te_uu = diag_remove(te_vortex_u).transpose()  # cener of vortices, remove diagonal and transpose
    te_u = np.tile(te_vortex_u, (len(te_vortex_u) - 1, 1))
    d3 = sum(-1j * te_ss / (2 * np.pi) * (1 / (te_u - te_uu)))

    te_ss = np.tile(te_vortex_strength, (len(te_vortex_strength), 1)).transpose()
    te_uu = np.conjugate(np.tile(te_vortex_u, (len(te_vortex_u), 1)).transpose())
    te_u = np.tile(te_vortex_u, (len(te_vortex_u), 1))
    d4 = sum(-1j * te_ss / (2 * np.pi) * (1 / (te_u * (te_u - te_uu))))

    p = d1 + d2 + d3 + d4 + d4L

    vel_conj = p * dudz_te - 1j * te_vortex_strength * d2udz2_te / dudz_te / (4 * np.pi)
    te_vortex_vel = np.conj(vel_conj)
    te_vortex_z = te_vortex_z + te_vortex_vel * time_step

    te_vortex_v = [newton(te_vortex_v[index], 1e-8, 250, Gkn, radius, center_circle, te_vortex_z[index])
                   for index in range(len(te_vortex_z))]
    te_vortex_v = np.array(te_vortex_v)
    te_vortex_u = (te_vortex_v - center_circle) / radius

    # - move vortices - leading edge emmits
    if len(le_vortex_z) > 0:
        # velocity
        d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / le_vortex_u ** 2)
        # circulation
        d2 = -1j * circulation / (2 * np.pi * le_vortex_u)
        # trailing edge emits vortices
        te_ss = np.tile(te_vortex_strength, (len(te_vortex_strength), 1)).transpose()
        te_uu = np.conjugate(np.tile(te_vortex_u, (len(te_vortex_u), 1)).transpose())
        te_u = np.tile(te_vortex_u, (len(te_vortex_u), 1))
        d4T = sum(-1j * te_ss / (2 * np.pi) * (1 / (te_u * (te_u - te_uu)))) + \
              sum(-1j * te_ss / (2 * np.pi) * (1 / (te_u - te_uu)))

        # leading edge emits vortices
        le_ss = diag_remove(le_vortex_strength).transpose()  # strength of vortices, remove diagonal and transpose
        le_uu = diag_remove(le_vortex_u).transpose()  # cener of vortices, remove diagonal and transpose
        le_u = np.tile(le_vortex_u, (len(le_vortex_u) - 1, 1))
        d3 = sum(-1j * le_ss / (2 * np.pi) * (1 / (le_u - le_uu)))

        le_ss = np.tile(le_vortex_strength, (len(le_vortex_strength), 1)).transpose()
        le_uu = np.conjugate(np.tile(le_vortex_u, (len(le_vortex_u), 1)).transpose())
        le_u = np.tile(le_vortex_u, (len(le_vortex_u), 1))
        d4 = sum(-1j * le_ss / (2 * np.pi) * (1 / (le_u * (le_u - le_uu))))

        p = d1 + d2 + d3 + d4 + d4T

        vel_conj = p * dudz_le - 1j * le_vortex_strength * d2udz2_le / dudz_le / (4 * np.pi)
        le_vortex_vel = np.conj(vel_conj)
        le_vortex_z = le_vortex_z + le_vortex_vel * time_step

        le_vortex_v = [newton(le_vortex_v[index], 1e-8, 250, Gkn, radius, center_circle, le_vortex_z[index])
                       for index in range(len(le_vortex_z))]
        le_vortex_v = np.array(le_vortex_v)
        le_vortex_u = (le_vortex_v - center_circle) / radius

    current_time += time_step

write_array(circulation_list, te_vortex_strength, iterate_time_step, heading_file)
print('total time ', time.time() - start)
