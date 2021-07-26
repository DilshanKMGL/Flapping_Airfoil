# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 18:02:32 2020

@author: dilsh
"""

"""
np.sum(a,b,axis=p) p=0 - column, p=1 - row
"""
import numpy as np
import time
from _datetime import datetime as dt
import os


def multiply_dif_lists(list_1, list_2):
    """
    

    Parameters
    ----------
    list_1 : nd array
        1st list
    list_2 : nd array
        2nd list
    
    Returns
    -------
    nd array
        selected 2D array

    ex:-
    list_1 = np.array([1,2,3,4,5])
    list_2 = np.array([10,20,30])

    return
    list_3 =
    [[ 10  20  30]
    [ 20  40  60]
    [ 30  60  90]
    [ 40  80 120]
    [ 50 100 150]]

    """
    list_1 = np.repeat(list_1, len(list_2)).reshape((len(list_1), len(list_2)))
    list_1 = list_1 * list_2
    return list_1


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


def newton(x0, epsilon, max_iter, Gkn, radius, center_circle, equal_val):
    """
    :param x0: required point
    :param epsilon: required tolerance
    :param max_iter: iteration per calculation
    :param Gkn: polynomial coefficients
    :param radius: radius of the circle in v plane
    :param center_circle: center of the circle in v plane
    :param equal_val: guess value in v plane
    :return: corresponding point in v plane
    """
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


def get_trailing_edge_u(trailing_edge_v, center_circle, radius):
    """
    Parameters
    ----------
    trailing_edge_v : complex
        trailing edge position in v plane
    center_circle : complex
        center of the circle in v plane
    radius : float
        radius of the circle in u plane

    Returns
    -------
    complex
        trailing edge position in u plane
    """
    return complex((trailing_edge_v - center_circle) / radius)


def get_steady_state_circulation(velocity, aoa, radius, trailing_edge_u):
    """
    Parameters
    ----------
    velocity : float
        velocity value of free stream
    aoa : float
        angle of attack of free stram
    radius : float
        radius of the circle in v plane
    trailing_edge_u : complex
        trailing edge position in u plane

    Returns
    -------
    float
        steady state circulation value

    """
    u1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / trailing_edge_u ** 2)
    u2 = - 1j / trailing_edge_u / (2 * np.pi)

    steady_circulation = complex(- u1 / u2)
    return steady_circulation.real


def get_steady_state_lift(density, velocity, circulation):
    return -density * velocity * circulation


def get_steady_state_cl(circulation, velocity):
    return -2 * circulation / velocity


def get_plunging_velocity(pl_amp, pl_freq, c_time):
    pl_vel = -2 * 1j * pl_amp * np.pi * pl_freq * np.sin(2 * np.pi * pl_freq * c_time)
    return pl_vel


def get_free_vel(re_num, viscosity, density):
    return re_num * viscosity / density


def get_relative_vel(free_vel, free_aoa, pl_vel):
    rel_vel = np.abs(free_vel + pl_vel)
    rel_aoa = free_aoa + np.angle(free_vel + pl_vel)
    return rel_vel, rel_aoa


def get_sr_num(pl_amplitude, pl_frequency, free_velocity):
    return 2 * np.pi * pl_amplitude * pl_frequency / free_velocity


def get_new_vor_pos(trailing_edge_z, distance, angle):  # under develop
    new_vortex_position_z = complex(trailing_edge_z + distance * pow(np.e, -1j * angle))


def cal_circulation(radius, velocity, aoa, trailing_edge_u, te_vortex_strength, new_vortex_position_u, te_vortex_u):
    """

    :param radius: radius of the circle in v plane
    :param velocity: freestream velocity
    :param aoa: angle of attack in radians
    :param trailing_edge_u: corresponding trailing edge of airfoil in u plane
    :param te_vortex_strength: nd array of trailing edge shed vortices
    :param new_vortex_position_u: new position of the vortex
    :param te_vortex_u: nd array vortices position in u plane
    :return:
    circulation: value of the circulation
    s: sum of vortices strength
    """
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


def cal_force(iterate, velocity, aoa, aoa_local, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre, circulation, te_vortex_u,
              te_vortex_z, te_vortex_strength, time_step, density, radius, center_circle, Gkn):
    """

    :param aerofoil:
    :param aoa_local:
    :param iterate: current iteration value
    :param velocity: freestream value relative to the aerofoil
    :param aoa: freestream angle of attack relative to the aerofoil
    :param Iwx_pre: wake impulse x direction in previous iteration
    :param Iwy_pre: wake impulse y direction in previous iteration
    :param Ibvx_pre: bound impulse in x direction in previous iteration
    :param Ibvy_pre: bound impulse in y direction in previous iteration
    :param circulation: current circulaion value
    :param te_vortex_u: nd array vortices positions in u plane
    :param te_vortex_z: nd array vortices positions in z plane
    :param te_vortex_strength: nd array strength of vortices
    :param time_step: value of time step
    :param radius: radius of the circle in v plane
    :param center_circle: center of the circle in v plane
    :param density: density of the flow
    :return:
    Fwx: wake force in x direction
    Fwy: wake force in y direction
    Fbvx: bound force in x direction
    Fbvy: bound force in y direction
    force_x: force in x direction
    force_y: force in y direction
    lift: lift value
    drag: drag value
    cl: coefficient of lift
    cd: coefficient of drag
    Iwx_pre: wake impluse in x direction for next iteration
    Iwy_pre: wake impluse in y direction for next iteration
    Ibvx_pre: bound impulse in x direction for next iteration
    Ibvy_pre: bound impulse in y direction for next iteration
    """
    # N, radius, center_circle, trailing_edge_z, trailing_edge_v, Gkn, z_plane, v_plane, u_plane = read_data(aerofoil)
    # calculate wake vorticitym
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
    lift = force_y * np.cos(aoa_local) - force_x * np.sin(aoa_local)
    drag = force_x * np.cos(aoa_local) + force_y * np.sin(aoa_local)

    denom = density * np.power(velocity, 2)
    cl = 2 * lift / denom
    cd = 2 * drag / denom

    return cl, cd, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre


def cal_velocity(velocity, aoa, circulation, te_vor_u, te_vor_strength, Gkn, radius, req_point_u):
    power = np.arange(1, len(Gkn) + 1)
    dzdu = radius - sum(Gkn * power / req_point_u ** (power + 1))

    # free stream
    d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / req_point_u ** 2)
    # circulation
    d2 = -1j * circulation / (2 * np.pi * req_point_u)

    p1 = 1.0 / (req_point_u - te_vor_u)
    p2 = 1.0 / (req_point_u * (1.0 - req_point_u * np.conjugate(te_vor_u)))
    d4 = -1j * te_vor_strength / (2 * np.pi) * (p1 + p2)
    d4 = sum(d4)

    dwdu = d1 + d2 + d4
    velocity_point = np.conj(dwdu / dzdu)
    return velocity_point


def move_vortices(te_vortex_u, te_vortex_v, te_vortex_z, te_vortex_strength, velocity, aoa, circulation, aerofoil,
                  time_step):
    """

    :param aerofoil: aerofoil name
    :param te_vortex_u: nd array of vortices positions in u plane
    :param te_vortex_v: nd array of vortices positions in v plane
    :param te_vortex_z: nd array of vortices positions in z plane
    :param velocity: freestream velocity
    :param aoa: freestream angle of attack relative to aerofoil
    :param circulation: current circulation value
    :param te_vortex_strength: nd array strength of trailing edge vortices
    :param time_step: time step value
    :return:
    te_vortex_u: new vortex points in u plane
    te_vortex_v: new vortex points in v plane
    te_vortex_z: new vortex points in z plane
    """

    N, radius, center_circle, trailing_edge_z, trailing_edge_v, Gkn, z_plane, v_plane, u_plane = read_data(aerofoil)
    # - calculate derivative
    power = np.arange(1, N + 1)
    Gkn_coeff = np.tile(Gkn, (len(te_vortex_u), 1)).transpose()
    power_coeff = np.tile(power, (len(te_vortex_u), 1)).transpose()
    te_coeff = np.tile(te_vortex_u, (N, 1))

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
        # for real part
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


def read_data(aerofoil):
    heading = 'Aerofoil_data/' + str(aerofoil) + '_data.txt'
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
    file1.close()

    return N, r, center_circle, trailing_edge_z, trailing_edge_v, Gkn, z_plane, v_plane, u_plane


def make_position_file(airfoil, re_num, sr_number, density, viscosity, free_velocity, free_aoa, pl_amplitude,
                       pl_frequency, time_step, iteration, heading):
    file1 = open(heading, 'w')

    file1.write('airfoil\n' + str(airfoil) + '\n')
    file1.write('re_num\n' + str(re_num) + '\n')
    file1.write('strauhl number\n' + str(sr_number) + '\n')
    file1.write('density\n' + str(density) + '\n')
    file1.write('viscosity\n' + str(viscosity) + '\n')
    file1.write('free_velocity\n' + str(free_velocity) + '\n')
    file1.write('free_aoa\n' + str(free_aoa) + '\n')
    file1.write('pl_amplitude\n' + str(pl_amplitude) + '\n')
    file1.write('pl_frequency\n' + str(pl_frequency) + '\n')
    file1.write('time_step\n' + str(time_step) + '\n')
    file1.write('iteration\n' + str(iteration) + '\n')
    file1.write('te_vortex_z' + '\n')
    file1.write('First value is the iteration number' + '\n')
    file1.close()


def make_force_file(airfoil, re_num, sr_number, density, viscosity, free_velocity, free_aoa, pl_amplitude,
                    pl_frequency, time_step, iteration, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, 'w')

    file1.write('airfoil\n' + str(airfoil) + '\n')
    file1.write('re_num\n' + str(re_num) + '\n')
    file1.write('strauhl number\n' + str(sr_number) + '\n')
    file1.write('density\n' + str(density) + '\n')
    file1.write('viscosity\n' + str(viscosity) + '\n')
    file1.write('free_velocity\n' + str(free_velocity) + '\n')
    file1.write('free_aoa\n' + str(free_aoa) + '\n')
    file1.write('pl_amplitude\n' + str(pl_amplitude) + '\n')
    file1.write('pl_frequency\n' + str(pl_frequency) + '\n')
    file1.write('time_step\n' + str(time_step) + '\n')
    file1.write('iteration\n' + str(iteration) + '\n')
    file1.write('First value is the iteration number' + '\n')
    file1.write('Force values - Cl, Cd\n')
    file1.close()


def update_position_file(te_vortex_z, iteration, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write(str(iteration) + ' ')
    for value in te_vortex_z:
        file1.write(str(value) + ' ')
    file1.write('\n')
    file1.close()


def update_force_file(iterate, cl, cd, heading):
    file1 = open(heading, "a+")
    file1.write(str(iterate) + ' ')  # iteration number
    file1.write(str(cl) + ' ' + str(cd) + '\n')
    file1.close()


def write_strength(circ_list, te_vor_str, end_time, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write('circulation\n')
    for value in circ_list:
        file1.write(str(value) + ' ')
    file1.write('\n')
    file1.write('te_vortex_strength\n')
    for value in te_vor_str:
        file1.write(str(value) + ' ')
    file1.write('\n')
    file1.write('complete time\n')
    file1.write(str(end_time) + 's\n')
    file1.close()


# velocity definition has been changed
def main(sr_number, reduced_frequency, free_velocity, amplitude, time_step):
    start = time.time()
    print("This is the main function")

    # aerofoil
    aerofoil = 'NACA2412'

    # free stream parameters
    re_num = 1.0e6
    density = 1.225  # km/m^3
    viscosity = 1.789e-5  # Pa.s
    free_aoa_deg = 0.0  # degrees
    # velocity, aoa calculation
    # free_velocity = 10  # get_free_vel(re_num, viscosity, density)
    free_aoa = np.deg2rad(free_aoa_deg)

    # plunging parameters
    # sr_number = 0.48
    # reduced_frequency = 0.753982237
    pl_frequency = reduced_frequency * free_velocity / np.pi
    pl_amplitude = sr_number * free_velocity / pl_frequency
    # pl_frequency = 0.1909859317
    # pl_amplitude = 1.570796327
    # sr_number = get_sr_num(pl_amplitude, pl_frequency, free_velocity)
    # reduced_frequency = np.pi * pl_frequency / free_velocity

    # time step
    # time_step = 0.0005
    current_time = 0.0
    iteration = 2500  # use maximum iterations as 3500. otherwise it will get memory error

    # calculate new vortex position
    distance = 0.0005
    angle_deg = 0.0
    error_margine = 0.1  # percentage error

    # data_storage
    circulation_list = np.array([])
    te_vortex_strength = np.array([])
    te_vortex_z = np.array([])
    te_vortex_v = np.array([])
    te_vortex_u = np.array([])
    le_vortex_strength = np.array([])
    le_vortex_z = np.array([])
    le_vortex_v = np.array([])
    le_vortex_u = np.array([])
    cl_array = np.array([])
    cd_array = np.array([])
    plunging_dis_array = np.array([])
    iterate_time_step = np.array([])

    # force parameters
    Iwx_pre = 0
    Iwy_pre = 0
    Ibvx_pre = 0
    Ibvy_pre = 0

    # data read
    N, radius, center_circle, trailing_edge_z, trailing_edge_v, Gkn, z_plane, v_plane, u_plane = read_data(aerofoil)
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)
    # steady state calculation
    # trailing_edge_u = get_trailing_edge_u(trailing_edge_v, center_circle, radius)
    # steady_circulation = get_steady_state_circulation(free_velocity, free_aoa, radius, trailing_edge_u)
    # steady_cl = get_steady_state_cl(steady_circulation, free_velocity)
    # steady_lift = get_steady_state_lift(density, free_velocity, steady_circulation)

    # calculate vortex position
    angle = np.deg2rad(angle_deg)
    new_vor_pos_z = trailing_edge_z + distance * pow(np.e, -1j * angle)
    new_vor_pos_z = complex(new_vor_pos_z)

    # make file
    if not os.path.exists('Results'):
        os.mkdir('Results')
    ctime = dt.now().strftime("%Y-%m-%d %H.%M.%S")
    path_dir = 'Results/' + ctime + ' ' + aerofoil + ' plunging'
    path_dir = 'Results/' + aerofoil + '  TESM Re-' + str(re_num) + ' ,St-' + str(round(sr_number, 4)) + ' ,k-' + \
               str(round(reduced_frequency, 4)) + ' ,aoa-' + str(free_aoa_deg)
    os.mkdir(path_dir)

    position_file = path_dir + '/result_file.txt'
    is_position = False
    make_position_file(aerofoil, re_num, sr_number, density, viscosity, free_velocity, free_aoa_deg, pl_amplitude,
                       pl_frequency, time_step, iteration, position_file)

    force_file = path_dir + '/force_file.txt'
    make_force_file(aerofoil, re_num, sr_number, density, viscosity, free_velocity, free_aoa_deg, pl_amplitude,
                    pl_frequency, time_step, iteration, force_file)

    itr_time = time.time()
    for iterate in range(iteration):
        # calculate velocity
        velocity = free_velocity
        aoa = free_aoa
        aoa_local = aoa
        if pl_amplitude != 0.0 and pl_frequency != 0.0:
            pl_velocity = get_plunging_velocity(pl_amplitude, pl_frequency, current_time)
            velocity = np.abs(free_velocity + pl_velocity)
            aoa = free_aoa + np.angle(free_velocity + pl_velocity)
            aoa_local = -aoa

        cont_iteration = 1
        while True:
            search_point = center_circle + radius
            new_vor_pos_v = complex(newton(search_point, 1e-6, 100, Gkn, radius, center_circle, new_vor_pos_z))
            new_vor_pos_u = complex((new_vor_pos_v - center_circle) / radius)

            req_point_z = trailing_edge_z + (new_vor_pos_z - trailing_edge_z) / 2.0
            # req_point_z = (new_vor_pos_z + trailing_edge_z) / 2.0
            req_point_v = complex(newton(new_vor_pos_v, 1e-8, 100, Gkn, radius, center_circle, req_point_z))
            req_point_u = (req_point_v - center_circle) / radius

            # calculate circulation
            circulation, s = cal_circulation(radius, velocity, aoa, trailing_edge_u, te_vortex_strength, new_vor_pos_u,
                                             te_vortex_u)
            velocity_point = cal_velocity(velocity, aoa, circulation, te_vortex_u, te_vortex_strength, Gkn, radius,
                                          req_point_u)
            distance_point = complex(trailing_edge_z + velocity_point * time_step)
            error = abs((abs(distance_point) - abs(new_vor_pos_z)) / abs(distance_point) * 100)

            if error > error_margine:
                new_vor_pos_z = (new_vor_pos_z + distance_point) / 2.0
                cont_iteration += 1
            else:
                # print('iterate', iterate)
                # print('error', error)
                # print('error_margine', error_margine)
                # print('new_vor_pos_z', new_vor_pos_z)
                # print('req_point_z', req_point_z)
                # print('distance_point', distance_point)
                # print(iterate, 'error', error, 'cont_iteration', cont_iteration)
                break

        # print(circulation)
        # data update
        circulation_list = np.append(circulation_list, [circulation])
        te_vortex_strength = np.append(te_vortex_strength, [- s - circulation])
        te_vortex_z = np.append(te_vortex_z, [new_vor_pos_z])
        te_vortex_v = np.append(te_vortex_v, [new_vor_pos_v])
        te_vortex_u = np.append(te_vortex_u, [new_vor_pos_u])

        # leading edge
        le_z = min(z_plane)
        search_point = center_circle - radius
        le_v = complex(newton(search_point, 1e-6, 100, Gkn, radius, center_circle, le_z))
        le_u = complex((le_v - center_circle) / radius)

        vel_leading = cal_velocity(velocity, aoa, circulation, te_vortex_u, te_vortex_strength, Gkn, radius, le_u)

        # update position file
        if is_position:
            update_position_file(te_vortex_z, iterate, position_file)

        # calculate forces
        cl, cd, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre = \
            cal_force(iterate, velocity, aoa, aoa_local, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre, circulation, te_vortex_u,
                      te_vortex_z, te_vortex_strength, time_step, density, radius, center_circle, Gkn)

        cl_array = np.append(cl_array, [cl])
        cd_array = np.append(cd_array, [cd])

        # update force file
        update_force_file(iterate, cl, cd, force_file)

        # move vortices
        te_vortex_u, te_vortex_v, te_vortex_z = move_vortices(te_vortex_u, te_vortex_v, te_vortex_z, te_vortex_strength,
                                                              velocity, aoa, circulation, aerofoil, time_step)

        current_time += time_step
        if iterate % 100 == 0:
            print('time for', iterate, " iteration:", time.time() - itr_time, 's')
            itr_time = time.time()
    end_time = time.time() - start
    print('total time ', time.time() - start)
    write_strength(circulation_list, te_vortex_strength, end_time, position_file)


sr_number_list = [0.48, 0.36, 0.32, 0.3, 0.2, 0.16, 0.12, 0.49, 0.34, 0.17]
frequency_list = [2.4, 1.8, 1.6, 1.5, 1.0, 0.8, 0.6, 2.45, 1.7, 0.85]
free_velocity_outer = 10
amplitude_outer = 1.0 * 2

for index in range(len(sr_number_list)):
    reduced_frequency_outer = np.pi * frequency_list[index] / free_velocity_outer
    time_step_outer = 1 / frequency_list[index] / 2500.0 * 2.3

    if __name__ == "__main__":
        main(sr_number_list[index], reduced_frequency_outer, free_velocity_outer, amplitude_outer, time_step_outer)
