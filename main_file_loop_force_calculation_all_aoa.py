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


def make_force_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
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
    file1.write('Force values - x and y, Fvx, Fvy, Fwx, Fwy, Fbvx, Fbvy\n')
    file1.close()


def make_mis_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
                  time_step, current_time, iteration, distance, angle, heading, type_name):
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
    file1.write(type_name + '\n')
    file1.close()


def update_file(te_vortex_z, iteration, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write('te_vortex_z - iteration ' + str(iteration + 1) + '\n' + str(list(te_vortex_z)) + '\n')
    file1.close()


def update_force_file(force_x, force_y, Fvx, Fvy, Fwx, Fwy, Fbvx, Fbvy, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write(str(force_x) + ' ' + str(force_y) + ' ' + str(Fvx) + ' ' + str(Fvy) + ' ' + str(Fwx) + ' ' +
                str(Fwy) + ' ' + str(Fbvx) + ' ' + str(Fbvy) + '\n')
    file1.close()


def update_mis_file(para1, heading, value):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    if len(list(para1)) > 1:
        for index in range(len(para1)):
            if value == 0:
                # wrtite complex numbers directly
                file1.write(str(para1[index]) + ' ')
            elif value == 1:
                # write real and imaginary parts seperately
                file1.write(str(para1[index].real) + ' ' + str(para1[index].imag) + ' ')
        file1.write('\n')
    else:
        file1.write(str(para1) + '\n')
    file1.close()


def write_array(circulation, te_vortex_strength, iterate_time_step, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
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


def initialize_field():
    # ------ calculate trailing edge position
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)

    # ------ calculate plunging parameters
    plunging_dis = 1j * pl_amplitude * np.sin(2 * np.pi * pl_frequency * current_time)
    plunging_vel = 2 * 1j * pl_amplitude * np.pi * pl_frequency * np.cos(2 * np.pi * pl_frequency * current_time)

    # --- calculate velocity
    velocity = np.abs(free_velocity + plunging_vel)
    aoa = free_aoa + np.angle(free_velocity + plunging_vel)

    # --- calculate velocity
    # velocity = free_velocity
    # aoa = free_aoa

    # ------ calculate new vortex position
    new_vortex_position_z = trailing_edge_z + distance * pow(np.e, -1j * angle)
    new_vortex_position_z = complex(new_vortex_position_z)
    search_point = center_circle + radius
    new_vortex_position_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, new_vortex_position_z))
    new_vortex_position_u = complex((new_vortex_position_v - center_circle) / radius)
    return trailing_edge_u, velocity, aoa, new_vortex_position_z, new_vortex_position_v, new_vortex_position_u


def calculate_circulation():
    # ------ create function to calculate circulation
    s = sum(te_vortex_strength)
    # velocity calculation
    d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / trailing_edge_u ** 2)
    # circulation
    d2 = -1j / (2 * np.pi * trailing_edge_u)
    # newly sheded vortex
    p1 = 1.0 / (trailing_edge_u - new_vortex_position_u)
    p2 = 1.0 / (trailing_edge_u * (1.0 - trailing_edge_u * np.conj(new_vortex_position_u)))
    d3 = -1j / (2 * np.pi) * (p1 + p2)

    # previously shed vortices
    p1 = 1.0 / (trailing_edge_u - te_vortex_u)
    p2 = 1.0 / (trailing_edge_u * (1.0 - trailing_edge_u * np.conjugate(te_vortex_u)))
    d4 = -1j * te_vortex_strength / (2 * np.pi) * (p1 + p2)
    d4 = sum(d4)

    circulation = complex((s * d3 - d1 - d4) / (d2 - d3)).real

    return circulation, s


def move_vortices(iterate_time_step, te_vortex_u, te_vortex_v, te_vortex_z):
    # - calculate derivative
    power = np.arange(1, len(Gkn) + 1)
    Gkn_coeff = np.tile(Gkn, (len(te_vortex_u), 1)).transpose()
    power_coeff = np.tile(power, (len(te_vortex_u), 1)).transpose()
    te_coeff = np.tile(te_vortex_u, (len(Gkn), 1))

    dudv = 1 / radius
    dzdv = 1.0 - sum(Gkn_coeff * power_coeff / radius / te_coeff ** (power_coeff + 1))
    dvdz = 1 / dzdv
    dudz = dudv * dvdz

    d2zdv2 = sum(Gkn_coeff * power_coeff * (power_coeff + 1) / radius ** 2 / te_coeff ** (power_coeff + 2))
    d2udz2 = - dudv * d2zdv2 / dzdv ** 3

    # ------ update file
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
        te_ss = diag_remove(te_vortex_strength).transpose()  # strength of vortices, remove diagonal and transpose
        te_uu = diag_remove(te_vortex_u).transpose()  # cener of vortices, remove diagonal and transpose
        te_u = np.tile(te_vortex_u, (len(te_vortex_u) - 1, 1))
        d3 = sum(-1j * te_ss / (2 * np.pi) * (1 / (te_u - te_uu)))

        te_ss = np.tile(te_vortex_strength, (len(te_vortex_strength), 1)).transpose()
        te_uu = np.conjugate(np.tile(te_vortex_u, (len(te_vortex_u), 1)).transpose())
        te_u = np.tile(te_vortex_u, (len(te_vortex_u), 1))
        d4 = sum(-1j * te_ss / (2 * np.pi) * (1 / (te_u * (1.0 - te_u * te_uu))))

        p += d3 + d4

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

    return te_vortex_u, te_vortex_v, te_vortex_z


def calcualte_force(iterate, Gkn, velocity, aoa, Iwvx_pre, Iwvy_pre, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre, circulation,
                        te_vortex_u, te_vortex_z, te_vortex_strength):
    # calculate by velocity
    I0 = velocity * np.exp(1j * aoa) * time_step * iteration
    Iwvx = I0.real
    Iwvy = I0.imag

    # calculate wake vorticity
    Iwx = sum(te_vortex_z.imag * te_vortex_strength)
    Iwy = -sum(te_vortex_z.real * te_vortex_strength)
    # airfoil_distance = velocity * np.exp(1j * aoa) * time_step * iteration
    # Iwx = sum((te_vortex_z-airfoil_distance).imag * te_vortex_strength)
    # Iwy = -sum((te_vortex_z-airfoil_distance).real * te_vortex_strength)

    # calculate bound vorticity
    num_div = 100
    phi = 2 * np.pi / num_div  # divisable number
    angle = np.array([n * phi for n in range(num_div)])
    circle_point_u = np.exp(1j * angle)

    n = np.arange(1, len(Gkn) + 1)
    circle_point_u_1 = np.tile(circle_point_u, (len(Gkn), 1))
    Gkn_1 = np.tile(Gkn, (len(circle_point_u), 1)).transpose()
    n_1 = np.tile(n, (len(circle_point_u), 1)).transpose()

    airfoil_point = circle_point_u * radius + center_circle + sum(Gkn_1 / (circle_point_u_1 ** n_1))

    # - calculate derivative
    power = np.arange(1, len(Gkn) + 1)
    Gkn_coeff = np.tile(Gkn, (len(circle_point_u), 1)).transpose()
    power_coeff = np.tile(power, (len(circle_point_u), 1)).transpose()
    te_coeff = np.tile(circle_point_u, (len(Gkn), 1))

    dzdv = 1.0 - sum(Gkn_coeff * power_coeff / radius / te_coeff ** (power_coeff + 1))

    ''' calculate dwdv '''
    # velocity
    d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / circle_point_u ** 2)
    # circulation
    d2 = -1j * circulation / (2 * np.pi * circle_point_u)
    # vortices
    point1 = np.tile(circle_point_u, (len(te_vortex_strength), 1))
    point2 = np.tile(te_vortex_u, (len(circle_point_u), 1)).transpose()
    strength = np.tile(te_vortex_strength, (len(circle_point_u), 1)).transpose()

    p1 = 1.0 / (point1 - point2)
    p2 = 1.0 / (point1 * (1.0 - point1 * np.conjugate(point2)))
    d3 = -1j * strength / (2 * np.pi) * (p1 + p2)
    d3 = sum(d3)
    dwdv = d1 + d2 + d3

    # inner = np.imag((dwdv - velocity * np.exp(-1j * aoa) * dzdv) * np.exp(1j * angle)) * phi
    inner = np.imag(dwdv / dzdv * np.exp(1j * angle)) * phi
    Ibvx = - sum(inner * airfoil_point.imag)
    Ibvy = sum(inner * airfoil_point.real)

    '''mis file section'''
    if mis_file:
        # tangential_velocity = (dwdv - velocity * np.exp(-1j * aoa) * dzdv) * np.exp(1j * angle)
        # update_mis_file(tangential_velocity, heading_mis_file)
        velocity_air = dwdv/dzdv
        update_mis_file(velocity_air, heading_mis_file, 1)

    # force calculation
    if iterate != 0:
        Fvx = - (Iwvx - Iwvx_pre) / time_step * density
        Fvy = - (Iwvy - Iwvy_pre) / time_step * density
        Fwx = - (Iwx - Iwx_pre) / time_step * density
        Fwy = - (Iwy - Iwy_pre) / time_step * density
        Fbvx = - (Ibvx - Ibvx_pre) / time_step
        Fbvy = - (Ibvy - Ibvy_pre) / time_step
    else:
        Fvx = 0
        Fvy = 0
        Fwx = 0
        Fwy = 0
        Fbvx = 0
        Fbvy = 0

    Iwvx_pre = Iwvx
    Iwvy_pre = Iwvy
    Iwx_pre = Iwx
    Iwy_pre = Iwy
    Ibvx_pre = Ibvx
    Ibvy_pre = Ibvy

    return Fvx, Fvy, Fwx, Fwy, Fbvx, Fbvy, Iwvx_pre, Iwvy_pre, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre


aoa_list = np.arange(1.0, 20.0, 1)

for aoa_new in aoa_list:
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
    free_aoa = aoa_new
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
    current_time = 0.00
    iteration = 1000

    # ----- write in a file
    heading_file = 'Transient_solution_results/' + 'result_file_' + airfoil + '_aoa_' + str(aoa_new) + '.txt'
    make_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
              time_step, current_time, iteration, distance, angle, heading_file)
    heading_force_file = 'Transient_solution_results/' + 'force_file_' + airfoil + '_aoa_' + str(aoa_new) + '.txt'
    make_force_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
                    time_step, current_time, iteration, distance, angle, heading_force_file)
    type_name = 'velocity_on_the_airfoil_'
    # type_name = 'tangential_velocity_'
    # type_name = 'velocity_on_the_airfoil_values_'
    # type_name = 'tangential_velocity_values_'
    heading_mis_file = 'Transient_solution_results/' + 'mis_file_' + type_name + airfoil + '.txt'
    mis_file = False
    if mis_file:
        make_mis_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
                      time_step, current_time, iteration, distance, angle, heading_mis_file, type_name)
    print(airfoil)

    Iwvx_pre = 0
    Iwvy_pre = 0

    Iwx_pre = 0
    Iwy_pre = 0

    Ibvx_pre = 0
    Ibvy_pre = 0

    # ------ iteration code
    for iterate in range(iteration):
        if iterate % 100 == 0:
            print('Iteration - ' + str(iterate))

        trailing_edge_u, velocity, aoa, new_vortex_position_z, new_vortex_position_v, new_vortex_position_u = \
            initialize_field()
        circulation, s = calculate_circulation()

        circulation_list = np.append(circulation_list, [circulation])
        te_vortex_strength = np.append(te_vortex_strength, [- s - circulation])
        te_vortex_z = np.append(te_vortex_z, [new_vortex_position_z])
        te_vortex_v = np.append(te_vortex_v, [new_vortex_position_v])
        te_vortex_u = np.append(te_vortex_u, [new_vortex_position_u])

        # - calculate forces
        Fvx, Fvy, Fwx, Fwy, Fbvx, Fbvy, Iwvx_pre, Iwvy_pre, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre = \
            calcualte_force(iterate, Gkn, velocity, aoa, Iwvx_pre, Iwvy_pre, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre,
                            circulation, te_vortex_u, te_vortex_z, te_vortex_strength)

        force_x = Fwx + Fbvx + Fvx
        force_y = Fwy + Fbvy + Fvy
        update_force_file(force_x, force_y, Fvx, Fvy, Fwx, Fwy, Fbvx, Fbvy, heading_force_file)
        iterate_time_step = np.append(iterate_time_step, [time.time() - iterate_time])

        te_vortex_u, te_vortex_v, te_vortex_z = move_vortices(iterate_time_step, te_vortex_u, te_vortex_v, te_vortex_z)
        current_time += time_step

    write_array(circulation_list, te_vortex_strength, iterate_time_step, heading_file)
    print('total time ', time.time() - start)