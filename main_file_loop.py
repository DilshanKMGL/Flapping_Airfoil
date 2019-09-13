import numpy as np
import time


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
    file1.write('te_vortex_z - iteration ' + str(iteration+1) + '\n' + str(list(te_vortex_z)) + '\n')


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
free_velocity = 5
free_aoa = 0.0
free_aoa = np.deg2rad(free_aoa)
# ------ plunging parameters
pl_amplitude = 0
pl_frequency = 0
# ------ pitching parameters
pi_amplitude = 0
pi_frequency = 0
# ------ new vortex
distance = 0.005
angle = -5
angle = np.deg2rad(angle)
# ------ data store
circulation_list = np.array([])
te_vortex_strength = np.array([])
te_vortex_z = np.array([])
te_vortex_v = np.array([])
te_vortex_u = np.array([])
iterate_time_step = np.array([])
# ------ time step
time_step = 0.05
current_time = 0.00
iteration = 200

# ------ derivatives

# ----- write in a file
make_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
             time_step, current_time, iteration, distance, angle)

# ------ iteration code
for iterate in range(iteration):
    print('Iteration - ' + str(iterate + 1))

    # ------ calculate trailing edge position
    search_point = center_circle + radius
    trailing_edge_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, trailing_edge_z))
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)
    # --- calculate velocity
    velocity = free_velocity
    aoa = free_aoa

    # ------ calculate new vortex position
    new_vortex_position_z = trailing_edge_z + distance * pow(np.e, -1j * angle)
    new_vortex_position_z = complex(new_vortex_position_z)
    search_point = center_circle + radius
    new_vortex_position_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, new_vortex_position_z))
    new_vortex_position_u = complex((new_vortex_position_v - center_circle) / radius)
    te_vortex_z = np.append(te_vortex_z, [new_vortex_position_z])
    te_vortex_v = np.append(te_vortex_v, [new_vortex_position_v])
    te_vortex_u = np.append(te_vortex_u, [new_vortex_position_u])

    # ------ create function to calculate circulation
    s = sum(te_vortex_strength)
    d1 = -1j / (2 * np.pi * trailing_edge_u)
    d2 = 1j * ((1 / (trailing_edge_u - new_vortex_position_u)) +
               (1 / (trailing_edge_u * (1 - trailing_edge_u * new_vortex_position_u)))) / (2 * np.pi)
    d3 = velocity * ((1 / pow(np.e, 1j * aoa)) - (pow(np.e, 1j * aoa) / (trailing_edge_u ** 2)))

    d4 = - 1j * te_vortex_strength * ((1 / (trailing_edge_u - te_vortex_u[:-1])) +
                                      (1 / (trailing_edge_u * (1 - trailing_edge_u * te_vortex_u[:-1])))) / (2 * np.pi)
    d4 = sum(d4)
    circulation = complex((- (s * d2 + d3 + d4) / (d1 + d2))).real
    circulation_list = np.append(circulation_list, [circulation])
    te_vortex_strength = np.append(te_vortex_strength, [-s - circulation])

    # - calculate derivative
    # te_vortex_v = te_vortex_u * radius + center_circle
    svc = te_vortex_v - center_circle
    svc = np.array(list(svc) * len(Gkn))
    svc = svc.reshape(len(Gkn), len(te_vortex_u))
    svc = np.transpose(svc)

    Gkncoeff = np.array(list(Gkn) * len(te_vortex_v))
    Gkncoeff = Gkncoeff.reshape(len(te_vortex_u), len(Gkn))

    power = np.arange(1, len(Gkn) + 1)
    power = np.array(list(power) * len(te_vortex_u))
    power = power.reshape(len(te_vortex_u), len(Gkn))

    dzdv = 1 - Gkncoeff * power * radius ** power / svc ** (power + 1)
    d2zdv2 = power * (power + 1) * Gkncoeff * radius ** power / svc ** (power + 2)
    dvdu = radius
    dudz = 1 / (dvdu * dzdv)
    d2udz2 = - d2zdv2 / dzdv ** 3

    # ------ iterate time
    iterate_time_step = np.append(iterate_time_step, [time.time() - iterate_time])
    update_file(te_vortex_z, iterate)

    # ------ move vortices
    p1 = -1j * circulation / te_vortex_u / (2 * np.pi)
    p2 = velocity * ((1 / pow(np.e, 1j * aoa)) - (pow(np.e, 1j * aoa) / (te_vortex_u ** 2)))

    u = te_vortex_u.copy()
    u = np.array(list(u) * len(u))
    u = u.reshape(len(te_vortex_u), len(te_vortex_u))
    strided = np.lib.stride_tricks.as_strided
    s0, s1 = u.strides
    u = strided(u.ravel()[1:], shape=(len(u) - 1, len(u)), strides=(s0 + s1, s1)).reshape(len(u), -1)

    vc = te_vortex_u.copy()
    vc = np.array(list(vc) * (len(vc) - 1))
    vc = vc.reshape(len(te_vortex_u) - 1, len(te_vortex_u))
    vc = np.transpose(vc)

    vs = np.transpose(te_vortex_strength).copy()
    vs = np.array(list(vs) * (len(vs) - 1))
    vs = vs.reshape(len(te_vortex_strength) - 1, len(te_vortex_strength))
    vs = np.transpose(vs)

    p = p1+p2
    if len(u) > 1:
        p3 = - 1j * vs * ((1 / (u - vc)) + (1 / (u * (1 - u * vc)))) / (2 * np.pi)
        p3 = sum(np.transpose(p3))
        p += p3

    p = np.array(list(p) * len(Gkn))
    p = p.reshape(len(Gkn), len(te_vortex_u))
    p = np.transpose(p)

    te_vortex = np.array(list(te_vortex_u) * len(Gkn))
    te_vortex = te_vortex.reshape(len(Gkn), len(te_vortex_u))
    te_vortex = np.transpose(te_vortex)

    vel_conj = p * dudz - 1j * te_vortex * d2udz2 / dudz / (4 * np.pi)
    vel_conj = sum(np.transpose(vel_conj))
    te_vortex_vel = np.conj(vel_conj)
    te_vortex_z = te_vortex_z + te_vortex_vel * time_step

    te_vortex_v = [newton(te_vortex_v[index], 1e-8, 50, Gkn, radius, center_circle, te_vortex_z[index])
                   for index in range(len(te_vortex_z))]
    te_vortex_v = np.array(te_vortex_v)
    te_vortex_u = (te_vortex_v - center_circle) / radius

    current_time += time_step

write_array(circulation_list, te_vortex_strength, iterate_time_step)
print('total time ', time.time() - start)
