import numpy as np
import time
from _datetime import datetime as dt
import os

from matplotlib import pyplot as plt
import graph
import write_files


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


def initialize_field(center_circle, radius, free_velocity, free_aoa, pl_frequency, pl_amplitude, trailing_edge_z,
                     distance, angle, Gkn, trailing_edge_v, current_time, plunging_on):
    # ------ calculate trailing edge position
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)

    # ------ calculate plunging parameters
    plunging_vel = 0
    plunging_dis = 0
    if plunging_on:
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
    return trailing_edge_u, velocity, aoa, new_vortex_position_z, new_vortex_position_v, new_vortex_position_u, \
           plunging_dis, plunging_vel


def calculate_circulation(radius, velocity, aoa, trailing_edge_u, te_vortex_strength, new_vortex_position_u,
                          te_vortex_u):
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


def move_vortices(te_vortex_u, te_vortex_v, te_vortex_z, Gkn, center_circle, radius, velocity, aoa, circulation,
                  te_vortex_strength, time_step):
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


def calcualte_force(iterate, Gkn, velocity, aoa, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre, circulation,
                    te_vortex_u, te_vortex_z, te_vortex_strength, time_step, radius, center_circle, density, mis_file,
                    heading_mis_file):
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
        # update_mis_file(iterate, tangential_velocity, heading_mis_file)
        velocity_air = dwdv / dzdv
        write_files.update_mis_file(iterate, velocity_air, heading_mis_file, 1)

    # force calculation
    if iterate != 0:
        Fwx = - (Iwx - Iwx_pre) / time_step * density
        Fwy = - (Iwy - Iwy_pre) / time_step * density
        Fbvx = - (Ibvx - Ibvx_pre) / time_step
        Fbvy = - (Ibvy - Ibvy_pre) / time_step

    else:
        Fwx = 0
        Fwy = 0
        Fbvx = 0
        Fbvy = 0

    Iwx_pre = Iwx
    Iwy_pre = Iwy
    Ibvx_pre = Ibvx
    Ibvy_pre = Ibvy

    force_x = Fwx + Fbvx
    force_y = Fwy + Fbvy
    lift = force_y * np.cos(aoa) - force_x * np.sin(aoa)
    drag = force_x * np.cos(aoa) + force_y * np.sin(aoa)

    denom = density * np.power(velocity, 2)
    cl = 2 * lift / denom
    cd = 2 * drag / denom

    return Fwx, Fwy, Fbvx, Fbvy, force_x, force_y, lift, drag, cl, cd, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre


def steady_state_circulation(velocity, aoa, trailing_edge_v, center_circle, radius):
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)
    u1 = velocity * radius * \
         (np.exp(-1j * aoa) - np.exp(1j * aoa) / trailing_edge_u ** 2)
    u2 = - 1j / trailing_edge_u / (2 * np.pi)

    steady_circulation = complex(- u1 / u2)
    return steady_circulation.real


def main():
    start = time.time()
    iterate_time = start
    ctime = dt.now().strftime("%Y-%m-%d %H.%M.%S")
    # ------ airfoil data
    airfoil = 'NACA0012'
    N, radius, center_circle, trailing_edge_z, trailing_edge_v, Gkn, z_plane, v_plane, u_plane = read_data(airfoil)

    # ------ free stream velocity
    re_num = 7e5
    density = 1.225
    viscosity = 1.789e-5
    free_velocity = re_num * viscosity / density
    free_aoa = 4.0
    free_aoa = np.deg2rad(free_aoa)

    # ------ plunging parameters
    pl_amplitude = 1.0
    pl_frequency = 5.0
    strauhl_number = 2 * np.pi * pl_amplitude * pl_frequency / free_velocity
    print('strauhl number', strauhl_number)
    plunging_on = True
    if pl_frequency == 0.0 or pl_amplitude == 0.0:
        plunging_on = False
    # ------ pitching parameters
    pi_amplitude = 0
    pi_frequency = 0

    # ------ new vortex
    distance = 0.001
    angle = 0
    angle = np.deg2rad(angle)

    # ------ time step
    time_step = 0.005
    current_time = 0.00
    iteration = 100

    # ----- force file parameters
    Iwx_pre = 0
    Iwy_pre = 0
    Ibvx_pre = 0
    Ibvy_pre = 0

    # ------ data store
    circulation_list = np.array([])
    te_vortex_strength = np.array([])
    te_vortex_z = np.array([])
    te_vortex_v = np.array([])
    te_vortex_u = np.array([])
    cl_array = np.array([])
    cd_array = np.array([])
    plunging_dis_array = np.array([])
    iterate_time_step = np.array([])

    # ----- writing file activation and conditions
    main_file = True
    force_file = True
    mis_file = False
    xlwrite = False

    # plot_graph_condition = True
    time_delay = time_step / 100  # time delay between two graph update intervals
    plot_all_4 = False
    nrow = 0
    ncol = 0
    heading_list = []
    x_axis_title = []
    y_axis_title = []
    if plot_all_4:
        nrow = 2
        ncol = 2

        heading_list.append('Iteration vs. Time')
        x_axis_title.append('Iteration')
        y_axis_title.append('Time (ms)')

        heading_list.append('Cl vs. Time')
        y_axis_title.append('Cl')
        x_axis_title.append('Time (s)')

        heading_list.append('Cd vs. Time')
        y_axis_title.append('Cd')
        x_axis_title.append('Time (s)')
    fig, axs = plt.subplots(ncols=ncol, nrows=nrow)

    vortex_movement = False
    cl_time_grapgh = False
    cd_time_grapgh = False
    plunging_dis_graph = True

    # ----- make directory
    if not os.path.exists('Results'):
        os.mkdir('Results')
    path_dir = 'Results/' + ctime
    os.mkdir(path_dir)

    # ----- write in a file
    heading_file = path_dir + '/result_file.txt'
    if main_file:
        write_files.make_file(airfoil, re_num, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude,
                              pi_frequency, time_step, current_time, iteration, distance, angle, heading_file)

    # heading_force_file = path_dir + '/force_file_' + airfoil + '.txt'
    heading_force_file = path_dir + '/force_file.txt'
    if force_file:
        write_files.make_force_file(airfoil, re_num, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude,
                                    pi_frequency, time_step, current_time, iteration, distance, angle,
                                    heading_force_file)

    type_name = 'velocity_on_the_airfoil_'
    # type_name = 'tangential_velocity_'
    # type_name = 'velocity_on_the_airfoil_values_'
    # type_name = 'tangential_velocity_values_'
    # heading_mis_file = path_dir + '/mis_file_' + type_name + airfoil + '.txt'
    heading_mis_file = path_dir + '/mis_file_' + type_name + '.txt'
    if mis_file:
        write_files.make_mis_file(airfoil, re_num, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude,
                                  pi_frequency, time_step, current_time, iteration, distance, angle, heading_mis_file,
                                  type_name)
    print(airfoil)

    # ------ steady state solution
    steady_circulation = steady_state_circulation(free_velocity, free_aoa, trailing_edge_v, center_circle, radius)
    steady_cl = -2 * steady_circulation / free_velocity
    steady_lift = -density * free_velocity * steady_circulation
    # ------ transcient solution
    for iterate in range(iteration):
        if iterate % 100 == 0:
            print('Iteration - ' + str(iterate))

        trailing_edge_u, velocity, aoa, new_vortex_position_z, new_vortex_position_v, new_vortex_position_u, \
        plunging_dis, plunging_vel = initialize_field(center_circle, radius, free_velocity, free_aoa, pl_frequency,
                                                      pl_amplitude, trailing_edge_z, distance, angle, Gkn,
                                                      trailing_edge_v, current_time, plunging_on)
        circulation, s = calculate_circulation(radius, velocity, aoa, trailing_edge_u, te_vortex_strength,
                                               new_vortex_position_u, te_vortex_u)

        circulation_list = np.append(circulation_list, [circulation])
        te_vortex_strength = np.append(te_vortex_strength, [- s - circulation])
        te_vortex_z = np.append(te_vortex_z, [new_vortex_position_z])
        te_vortex_v = np.append(te_vortex_v, [new_vortex_position_v])
        te_vortex_u = np.append(te_vortex_u, [new_vortex_position_u])

        # update main file
        if main_file:
            write_files.update_file(te_vortex_z, iterate, heading_file)

        # - calculate forces
        Fwx, Fwy, Fbvx, Fbvy, force_x, force_y, lift, drag, cl, cd, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre = \
            calcualte_force(iterate, Gkn, velocity, aoa, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre, circulation,
                            te_vortex_u, te_vortex_z, te_vortex_strength, time_step, radius, center_circle, density,
                            mis_file, heading_mis_file)

        cl_array = np.append(cl_array, [cl])
        cd_array = np.append(cd_array, [cd])
        if iterate == 1:
            cl_array[0] = cl_array[1]
            cd_array[0] = cd_array[1]

        if force_file:
            write_files.update_force_file(iterate, Fwx, Fwy, Fbvx, Fbvy, force_x, force_y, lift, drag, cl, cd,
                                          heading_force_file)
        iterate_time_step = np.append(iterate_time_step, [time.time() - iterate_time])

        # plot grapgh
        islast = False
        if iterate == iteration - 1:
            islast = True

        x_limit_low = -time_step * 2
        x_limit_high = iteration * time_step + time_step * 2
        iterate_array = np.arange(0, iterate + 1, 1)
        time_cal = time_step * iterate_array
        steady_value_list = steady_cl * np.ones(iterate + 1)
        if plot_all_4:
            x_data = []
            y_data = []

            x_data.append(iterate_array)
            y_data.append(iterate_time_step)
            x_data.append(time_cal)
            y_data.append(cl_array)
            x_data.append(time_cal)
            y_data.append(cd_array)

            graph.plot_graph_all(axs, heading_list, x_axis_title, y_axis_title, x_data, y_data, time_delay, islast,
                                 steady_value_list, x_limit_high, x_limit_low, path_dir, plunging_on)

        if cl_time_grapgh:
            graph.cl_grapgh_plot(time_cal, cl_array, time_delay, islast, x_limit_low, x_limit_high, steady_value_list,
                                 path_dir, plunging_on)
        if cd_time_grapgh:
            graph.cd_grapgh_plot(time_cal, cd_array, time_delay, islast, x_limit_low, x_limit_high, path_dir)
        if plunging_on and plunging_dis_graph:
            plunging_dis = 1j * pl_amplitude * np.sin(2 * np.pi * pl_frequency * current_time)
            plunging_dis_array = np.append(plunging_dis_array, [plunging_dis.imag])
            graph.plunging_distance_plot(time_cal, plunging_dis_array, time_delay, islast, x_limit_low, x_limit_high,
                                         path_dir)
        if vortex_movement:
            graph.plot_airfoil(z_plane, te_vortex_z, te_vortex_strength, free_aoa, time_delay, islast)
        te_vortex_u, te_vortex_v, te_vortex_z = move_vortices(te_vortex_u, te_vortex_v, te_vortex_z, Gkn, center_circle,
                                                              radius, velocity, aoa, circulation, te_vortex_strength,
                                                              time_step)
        current_time += time_step

    if main_file:
        if not plunging_on:
            write_files.write_steady_circulation(steady_circulation, steady_lift, steady_cl, heading_file)
        write_files.write_array(circulation_list, te_vortex_strength, iterate_time_step, heading_file)
    print('total time ', time.time() - start)

    if xlwrite:
        write_files.create_excel_file(path_dir)


main()
