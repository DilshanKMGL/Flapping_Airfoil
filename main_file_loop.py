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


def update_file(te_vortex_z, iteration, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write('te_vortex_z - iteration ' + str(iteration + 1) + '\n' + str(list(te_vortex_z)) + '\n')


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
    strided = np.lib.stride_tricks.as_strided

    p = np.tile(p, (len(p), 1))
    m = p.shape[0]
    s0, s1 = p.strides
    p = strided(p.ravel()[1:], shape=(m - 1, m), strides=(s0 + s1, s1)).reshape(m, -1).transpose()
    return p


start = time.time()
iterate_time = start
# ------ airfoil data
airfoil = 'NACA2412'
N, radius, center_circle, trailing_edge_z, trailing_edge_v, Gkn, z_plane, v_plane, u_plane = read_data(airfoil)
# ------ free stream velocity
re_num = 1e6
density = 1.225
viscosity = 1.789e-5
free_velocity = 20#re_num * viscosity / density
free_aoa = 0.0
free_aoa = np.deg2rad(free_aoa)
# ------ plunging parameters
pl_amplitude = 0
pl_frequency = 0
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
iterate_time_step = np.array([])
# ------ time step
time_step = 0.005
# " if the time step > 0.001, sudden variation of vortex position"
current_time = 0.00
iteration = 500

heading_file = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
# ----- write in a file
make_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
          time_step, current_time, iteration, distance, angle, heading_file)

# ------ iteration code
for iterate in range(iteration):
    print('Iteration - ' + str(iterate + 1))

    # ------ calculate trailing edge position
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)

    # --- calculate velocity
    velocity = free_velocity
    aoa = free_aoa
    # print(velocity)
    # print(aoa)

    # ------ calculate new vortex position
    new_vortex_position_z = trailing_edge_z + distance * pow(np.e, -1j * angle)
    new_vortex_position_z = complex(new_vortex_position_z)
    search_point = center_circle + radius
    new_vortex_position_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, new_vortex_position_z))
    new_vortex_position_u = complex((new_vortex_position_v - center_circle) / radius)

    # ------ create function to calculate circulation
    s = sum(te_vortex_strength)
    # velocity calculation
    d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / trailing_edge_u ** 2)
    # circulation
    d2 = -1j / (2 * np.pi * trailing_edge_u)
    # newly sheded vortex
    power = np.arange(1, len(Gkn) + 1)
    p1 = (trailing_edge_u * radius + center_circle) + \
         sum(Gkn / trailing_edge_u ** power) - \
         new_vortex_position_u
    p2 = radius / trailing_edge_u + \
         sum(Gkn * (radius * trailing_edge_u / (radius - trailing_edge_u * center_circle)) ** power) - \
         new_vortex_position_u
    p3 = radius - sum(Gkn * power / trailing_edge_u ** (power + 1))
    p4 = - radius / trailing_edge_u ** 2 + \
         sum(Gkn * power * trailing_edge_u ** (power - 1) *
             (radius / (radius - trailing_edge_u * center_circle)) ** (power + 1))
    d3 = -1j / (2 * np.pi) * (p4 / p2 - p3 / p1)
    # previously shed vortices
    p1 = (trailing_edge_u * radius + center_circle) + \
         sum(Gkn / trailing_edge_u ** power) - \
         te_vortex_u
    p2 = radius / trailing_edge_u + \
         sum(Gkn * (radius * trailing_edge_u / (radius - trailing_edge_u * center_circle)) ** power) - \
         te_vortex_u
    d4 = 1j * te_vortex_strength / (2 * np.pi) * (p4 / p2 - p3 / p1)
    d4 = sum(d4)

    circulation = complex(-(d1 + s * d3 + d4) / (d2 + d3)).real

    circulation_list = np.append(circulation_list, [circulation])
    te_vortex_strength = np.append(te_vortex_strength, [- s - circulation])
    te_vortex_z = np.append(te_vortex_z, [new_vortex_position_z])
    te_vortex_v = np.append(te_vortex_v, [new_vortex_position_v])
    te_vortex_u = np.append(te_vortex_u, [new_vortex_position_u])

    # ----------------------------- #
    #       should be revised       #
    # ----------------------------- #
    # - calculate derivative

    Gkn_coeff = np.tile(Gkn, (len(te_vortex_u), 1)).transpose()
    power_coeff = np.tile(power, (len(te_vortex_u), 1)).transpose()
    te_u = np.tile(te_vortex_u, (len(Gkn), 1))

    dudv = 1 / radius
    dzdv = sum(1 - Gkn_coeff * power_coeff / radius / te_u ** (power_coeff + 1))
    dvdz = 1 / dzdv
    dudz = dudv * dvdz

    d2zdv2 = sum(Gkn_coeff * power_coeff * (power_coeff + 1) / radius ** 2 / te_u ** (power_coeff + 2))
    d2udz2 = - dudv * d2zdv2 / dzdv ** 3

    # ------ iterate time
    iterate_time_step = np.append(iterate_time_step, [time.time() - iterate_time])
    update_file(te_vortex_z, iterate, heading_file)

    # ------------------------------------- #
    #           move vortices               #
    # ------------------------------------- #
    # velocity
    d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / te_vortex_u ** 2)
    # circulation
    d2 = -1j * circulation / (2 * np.pi * te_vortex_u)

    p = d1 + d2
    if len(te_vortex_u) > 1:
        # complex potential for shed vortices
        te_uu = diag_remove(te_vortex_u)  # square matrix of te_vortex_u, remove diagonal and transpose

        temp_1 = (te_vortex_u * radius + center_circle) + sum(Gkn_coeff / te_u ** power_coeff)
        temp_1 = np.tile(temp_1, (len(te_vortex_u) - 1, 1))
        p1 = sum(temp_1 - te_uu)

        temp_1 = (radius / te_vortex_u) + \
                 sum(Gkn_coeff * (te_u * radius / (radius - te_u * center_circle)) ** power_coeff)
        temp_1 = np.tile(temp_1, (len(te_vortex_u) - 1, 1))
        p2 = sum(temp_1 - te_uu)

        p3 = radius - sum(power_coeff * Gkn_coeff / te_u ** (power_coeff + 1))
        p4 = -radius / te_vortex_u ** 2 + \
             sum(Gkn_coeff * power_coeff * te_vortex_u ** (power_coeff - 1) *
                 (radius / (radius - te_u * center_circle)) ** (power_coeff + 1))
        d3 = 1j * te_vortex_strength / (2 * np.pi) * (p4 / p2 - p3 / p1)

        p += d3

    # ------------------------------------------------- #
    #            move vortices                          #
    # ------------------------------------------------- #

    vel_conj = p * dudz - 1j * te_vortex_strength * d2udz2 / dudz / (4 * np.pi)
    te_vortex_vel = np.conj(vel_conj)
    te_vortex_z = te_vortex_z + te_vortex_vel * time_step

    te_vortex_v = [newton(te_vortex_v[index], 1e-8, 250, Gkn, radius, center_circle, te_vortex_z[index])
                   for index in range(len(te_vortex_z))]
    te_vortex_v = np.array(te_vortex_v)
    te_vortex_u = (te_vortex_v - center_circle) / radius

    current_time += time_step

write_array(circulation_list, te_vortex_strength, iterate_time_step, heading_file)
print('total time ', time.time() - start)
