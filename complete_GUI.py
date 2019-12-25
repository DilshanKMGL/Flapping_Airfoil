import numpy as np
import time
from datetime import datetime as dt
import os
# import xlsxwriter as xl
from tkinter import *
from tkinter.ttk import *

from matplotlib import pyplot as plt


# import graph
# import write_files


def make_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
              time_step, current_time, iteration, distance, angle, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, 'w')

    file1.write('airfoil\n' + str(airfoil) + '\n')
    file1.write('free_velocity\n' + str(free_velocity) + '\n')
    file1.write('free_aoa\n' + str(np.rad2deg(free_aoa)) + '\n')
    file1.write('pl_amplitude\n' + str(pl_amplitude) + '\n')
    file1.write('pl_frequency\n' + str(pl_frequency) + '\n')
    file1.write('pi_amplitude\n' + str(pi_amplitude) + '\n')
    file1.write('pi_frequency\n' + str(pi_frequency) + '\n')
    file1.write('time_step\n' + str(time_step) + '\n')
    file1.write('current_time\n' + str(current_time) + '\n')
    file1.write('iteration\n' + str(iteration) + '\n')
    file1.write('distance\n' + str(distance) + '\n')
    file1.write('angle\n' + str(angle) + '\n')
    file1.write('te_vortex_z' + '\n')
    file1.write('First value is the iteration number' + '\n')
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
    file1.write('First value is the iteration number' + '\n')
    file1.write('Force values - force_x, force_y, Fwx, Fwy, Fbvx, Fbvy, lift, drag, Cl, Cd\n')
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
    file1.write('First value is the iteration number' + '\n')
    file1.write(type_name[:-1] + '\n')
    file1.close()


def update_file(te_vortex_z, iteration, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write(str(iteration) + ' ')
    for value in te_vortex_z:
        file1.write(str(value) + ' ')
    file1.write('\n')
    file1.close()


def update_force_file(iterate, Fwx, Fwy, Fbvx, Fbvy, force_x, force_y, lift, drag, cl, cd, heading):
    file1 = open(heading, "a+")
    file1.write(str(iterate) + ' ')  # iteration number
    file1.write(str(force_x) + ' ' + str(force_y) + ' ' + str(Fwx) + ' ' + str(Fwy) + ' ' + str(Fbvx) + ' ' +
                str(Fbvy) + ' ' + str(lift) + ' ' + str(drag) + ' ' + str(cl) + ' ' + str(cd) + '\n')
    file1.close()


def update_mis_file(iteration, para1, heading, value):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write(str(iteration) + ' ')
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


def write_steady_circulation(steady_circulation, steady_lift, steady_cl, heading):
    file1 = open(heading, "a+")
    file1.write('steasy state circulation\n' + str(steady_circulation) + '\n')
    file1.write('steasy state lift\n' + str(steady_lift) + '\n')
    file1.write('steasy state cl\n' + str(steady_cl) + '\n')


def write_array(circulation, te_vortex_strength, iterate_time_step, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, "a+")
    file1.write('circulation\n')
    for value in circulation:
        file1.write(str(value) + ' ')
    file1.write('\n')
    file1.write('te_vortex_strength\n')
    for value in te_vortex_strength:
        file1.write(str(value) + ' ')
    file1.write('\n')
    file1.write('iteration step time\n')
    for value in iterate_time_step:
        file1.write(str(value) + ' ')
    file1.write('\n')
    file1.close()


# def create_excel_file(path_dir):
#     force_file = open(path_dir + '/force_file.txt')
#     result_file = open(path_dir + '/result_file.txt')
#     force_line = force_file.readlines()
#     result_line = result_file.readlines()
#     airfoil_name = force_line[1][:-1]
#     workbook = xl.Workbook(path_dir + '/Data.xlsx')
#     worksheet = workbook.add_worksheet(airfoil_name)
#     workbook.close()


def plot_graph_all(axs, heading_list, x_axis_title, y_axis_title, x_data, y_data, pause_time, islast, steady_state_cl,
                   x_limit_high, x_limit_low, path_to_save, plunging_on):
    # iteration vs time
    axs[0, 0].set_xlabel(x_axis_title[0])
    axs[0, 0].set_ylabel(y_axis_title[0])
    axs[0, 0].set_title(heading_list[0])
    axs[0, 0].plot(x_data[0], y_data[0], color='B', label='plot 1', linewidth=1)
    axs[0, 0].grid(True)
    # cl/cd
    axs[0, 1].set_xlabel(x_axis_title[1])
    axs[0, 1].set_ylabel(y_axis_title[1])
    axs[0, 1].set_title('Cl &Cd')
    axs[0, 1].set_xlim([x_limit_low, x_limit_high])
    axs[0, 1].plot(x_data[1], y_data[1], color='R', label='plot 1', linewidth=1)
    axs[0, 1].plot(x_data[2], y_data[2], color='B', label='plot 1', linewidth=1)
    if not plunging_on:
        axs[0, 1].plot(x_data[1], steady_state_cl, color='y', label='plot 1', linewidth=1)
    axs[0, 1].grid(True)
    # cl
    axs[1, 0].set_xlabel(x_axis_title[1])
    axs[1, 0].set_ylabel(y_axis_title[1])
    axs[1, 0].set_title(heading_list[1])
    axs[1, 0].set_xlim([x_limit_low, x_limit_high])
    axs[1, 0].plot(x_data[1], y_data[1], color='R', label='plot 1', linewidth=1)
    if not plunging_on:
        axs[1, 0].plot(x_data[1], steady_state_cl, color='y', label='plot 1', linewidth=1)
    axs[1, 0].grid(True)
    # cd
    axs[1, 1].set_xlabel(x_axis_title[2])
    axs[1, 1].set_ylabel(y_axis_title[2])
    axs[1, 1].set_title(heading_list[2])
    axs[1, 1].set_xlim([x_limit_low, x_limit_high])
    axs[1, 1].plot(x_data[2], y_data[2], color='B', label='Cl', linewidth=1)
    axs[1, 1].grid(True)

    plt.tight_layout()
    plt.pause(pause_time)

    if islast:
        plt.savefig(path_to_save + '/Figure.png', dpi=1280)
        plt.show()


def cl_grapgh_plot(x1, y1, pause_time, islast, x_limit_low, x_limit_high, steady, path_to_save, plunging_on):
    plt.figure(1)
    plt.clf()
    plt.xlabel('Time (s)')
    plt.ylabel('Cl')
    plt.title('Cl vs time')
    plt.xlim(x_limit_low, x_limit_high)
    plt.plot(x1, y1, color='r', label='transient cl', linewidth=1)
    if not plunging_on:
        plt.plot(x1, steady, color='b', label='steady state cl', linewidth=1)
    plt.grid(True)
    plt.tight_layout()
    plt.pause(pause_time)
    plt.legend()
    if islast:
        plt.savefig(path_to_save + '/Cl vs time.png')
        plt.show()


def cd_grapgh_plot(x1, y1, pause_time, islast, x_limit_low, x_limit_high, path_to_save):
    plt.figure(2)
    plt.clf()
    plt.xlabel('Time (s)')
    plt.ylabel('Cd')
    plt.title('Cd vs time')
    plt.xlim(x_limit_low, x_limit_high)
    plt.plot(x1, y1, color='r', linewidth=1)
    plt.grid(True, which='both')
    plt.tight_layout()
    plt.pause(pause_time)

    if islast:
        plt.savefig(path_to_save + '/Cd vs time.png')
        plt.show()


def plunging_distance_plot(x1, y1, pause_time, islast, x_limit_low, x_limit_high, path_to_save):
    plt.figure(3)
    plt.clf()
    plt.xlabel('Time (s)')
    plt.ylabel('distance')
    plt.title('Plunging wing movement')
    plt.xlim(x_limit_low, x_limit_high)
    plt.plot(x1, y1, color='r', linewidth=1)
    plt.grid(True, which='both')
    plt.tight_layout()
    plt.pause(pause_time)

    if islast:
        plt.savefig(path_to_save + '/plunging movement.png')
        plt.show()


def plot_airfoil(airfoil, vortex, strength, aoa, pause_time, islast):  # free stream initial aoa should there
    aoa = np.exp(-1j * aoa)
    vortex_pos = np.array([])
    vortex_neg = np.array([])
    for index in range(len(vortex)):
        if strength[index] < 0:
            vortex_neg = np.append(vortex_neg, [vortex[index]])
        else:
            vortex_pos = np.append(vortex_pos, [vortex[index]])
    vortex_pos = vortex_pos * aoa
    vortex_neg = vortex_neg * aoa
    airfoil = airfoil * aoa

    plt.figure(4)
    plt.clf()
    plt.axis('off')
    plt.title('vortex movement')
    plt.grid(False)
    plt.axis('equal')
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(vortex_pos.real, vortex_pos.imag, s=2, color='b')
    plt.scatter(vortex_neg.real, vortex_neg.imag, s=2, color='r')
    plt.plot(airfoil.real, airfoil.imag)

    plt.tight_layout()
    plt.pause(pause_time)
    if islast:
        # plt.savefig(path_to_save + '/plunging movement.png')
        plt.show()


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
    if plunging_on:
        # plunging_dis = 1j * pl_amplitude * np.sin(2 * np.pi * pl_frequency * current_time)
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
        update_mis_file(iterate, velocity_air, heading_mis_file, 1)

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


def main(airfoil1, re_num1, density1, viscosity1, free_aoa1, pl_amplitude1, pl_frequency1, time_step1, iteration1,
         new_vor_distance1, new_vor_angle1, result_file1, force_file1, excel_file1, grph_vortex_move1, grph_clt1,
         grph_cdt1, grph_pld1):
    start = time.time()
    iterate_time = start
    ctime = dt.now().strftime("%Y-%m-%d %H.%M.%S")
    # ------ airfoil data
    airfoil = airfoil1
    N, radius, center_circle, trailing_edge_z, trailing_edge_v, Gkn, z_plane, v_plane, u_plane = read_data(airfoil)

    # ------ free stream velocity
    re_num = re_num1
    density = density1
    viscosity = viscosity1
    free_velocity = re_num * viscosity / density
    free_aoa = free_aoa1
    free_aoa = np.deg2rad(free_aoa)

    # ------ plunging parameters

    pl_amplitude = pl_amplitude1
    pl_frequency = pl_frequency1
    strauhl_number = 2 * np.pi * pl_amplitude * pl_frequency / free_velocity
    print('strauhl number', strauhl_number)
    plunging_on = True
    if pl_amplitude == 0 or pl_frequency == 0:
        plunging_on = False

    # ------ pitching parameters
    pi_amplitude = 0
    pi_frequency = 0

    # ------ new vortex
    distance = new_vor_distance1
    angle = new_vor_angle1
    angle = np.deg2rad(angle)

    # ------ time step
    time_step = time_step1
    current_time = 0.00
    iteration = iteration1

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
    main_file = result_file1
    force_file = force_file1
    mis_file = False
    xlwrite = excel_file1

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

    vortex_movement = grph_vortex_move1
    cl_time_grapgh = grph_clt1
    cd_time_grapgh = grph_cdt1
    plunging_dis_graph = grph_pld1

    # ----- make directory
    if not os.path.exists('Results'):
        os.mkdir('Results')
    path_dir = 'Results/' + ctime
    os.mkdir(path_dir)

    # ----- write in a file
    heading_file = path_dir + '/result_file.txt'
    if main_file:
        make_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
                  time_step, current_time, iteration, distance, angle, heading_file)

    # heading_force_file = path_dir + '/force_file_' + airfoil + '.txt'
    heading_force_file = path_dir + '/force_file.txt'
    if force_file:
        make_force_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude,
                        pi_frequency, time_step, current_time, iteration, distance, angle,
                        heading_force_file)

    type_name = 'velocity_on_the_airfoil_'
    # type_name = 'tangential_velocity_'
    # type_name = 'velocity_on_the_airfoil_values_'
    # type_name = 'tangential_velocity_values_'
    # heading_mis_file = path_dir + '/mis_file_' + type_name + airfoil + '.txt'
    heading_mis_file = path_dir + '/mis_file_' + type_name + '.txt'
    if mis_file:
        make_mis_file(airfoil, free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude,
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

        trailing_edge_u, velocity, aoa, new_vortex_position_z, new_vortex_position_v, new_vortex_position_u = \
            initialize_field(center_circle, radius, free_velocity, free_aoa, pl_frequency, pl_amplitude,
                             trailing_edge_z, distance, angle, Gkn, trailing_edge_v, current_time, plunging_on)
        circulation, s = calculate_circulation(radius, velocity, aoa, trailing_edge_u, te_vortex_strength,
                                               new_vortex_position_u, te_vortex_u)

        circulation_list = np.append(circulation_list, [circulation])
        te_vortex_strength = np.append(te_vortex_strength, [- s - circulation])
        te_vortex_z = np.append(te_vortex_z, [new_vortex_position_z])
        te_vortex_v = np.append(te_vortex_v, [new_vortex_position_v])
        te_vortex_u = np.append(te_vortex_u, [new_vortex_position_u])

        # update main file
        if main_file:
            update_file(te_vortex_z, iterate, heading_file)

        # - calculate forces
        Fwx, Fwy, Fbvx, Fbvy, force_x, force_y, lift, drag, cl, cd, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre = \
            calcualte_force(iterate, Gkn, velocity, aoa, Iwx_pre, Iwy_pre, Ibvx_pre, Ibvy_pre, circulation,
                            te_vortex_u, te_vortex_z, te_vortex_strength, time_step, radius, center_circle, density,
                            mis_file,
                            heading_mis_file)

        cl_array = np.append(cl_array, [cl])
        cd_array = np.append(cd_array, [cd])
        if iterate == 1:
            cl_array[0] = cl_array[1]
            cd_array[0] = cd_array[1]

        if force_file:
            update_force_file(iterate, Fwx, Fwy, Fbvx, Fbvy, force_x, force_y, lift, drag, cl, cd,
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

            plot_graph_all(axs, heading_list, x_axis_title, y_axis_title, x_data, y_data, time_delay, islast,
                           steady_value_list, x_limit_high, x_limit_low, path_dir, plunging_on)

        if cl_time_grapgh:
            cl_grapgh_plot(time_cal, cl_array, time_delay, islast, x_limit_low, x_limit_high, steady_value_list,
                           path_dir, plunging_on)
        if cd_time_grapgh:
            cd_grapgh_plot(time_cal, cd_array, time_delay, islast, x_limit_low, x_limit_high, path_dir)
        if plunging_on and plunging_dis_graph:
            plunging_dis = 1j * pl_amplitude * np.sin(2 * np.pi * pl_frequency * current_time)
            plunging_dis_array = np.append(plunging_dis_array, [plunging_dis.imag])
            plunging_distance_plot(time_cal, plunging_dis_array, time_delay, islast, x_limit_low, x_limit_high,
                                   path_dir)
        if vortex_movement:
            plot_airfoil(z_plane, te_vortex_z, te_vortex_strength, free_aoa, time_delay, islast)
        te_vortex_u, te_vortex_v, te_vortex_z = move_vortices(te_vortex_u, te_vortex_v, te_vortex_z, Gkn, center_circle,
                                                              radius, velocity, aoa, circulation, te_vortex_strength,
                                                              time_step)
        current_time += time_step

    if main_file:
        if not plunging_on:
            write_steady_circulation(steady_circulation, steady_lift, steady_cl, heading_file)
        write_array(circulation_list, te_vortex_strength, iterate_time_step, heading_file)
    print('total time ', time.time() - start)

    # if xlwrite:
    #     create_excel_file(path_dir)


def button_command():
    # solution.main()
    airfoil = combo.get()
    re_num = float(text4.get())
    density = float(text5.get())
    viscosity = float(text6.get())
    free_aoa = float(text7.get())
    pl_amplitude = float(text9.get())
    pl_frequency = float(text10.get())
    time_step = float(text12.get())
    iteration = int(text13.get())
    new_vor_distance = float(text15.get())
    new_vor_angle = float(text16.get())

    result_file = float(chk_state1.get())
    force_file = float(chk_state2.get())
    excel_file = False  # float(chk_state3.get())

    grph_vortex_move = float(chk_state4.get())
    grph_clt = float(chk_state5.get())
    grph_cdt = float(chk_state6.get())
    grph_pld = float(chk_state7.get())

    main(airfoil, re_num, density, viscosity, free_aoa, pl_amplitude, pl_frequency, time_step, iteration,
         new_vor_distance, new_vor_angle, result_file, force_file, excel_file, grph_vortex_move, grph_clt,
         grph_cdt, grph_pld)


def exit_command():
    exit()


label_width = 25
text_width = 15
window = Tk()
window.title("DVM Application")

# column 1
label0 = Label(window, text=' ')
label0.grid(row=0, column=0)

label1 = Label(window, text='Aerofoil selection', font=('Arial Black', 10), width=label_width, anchor=W)
label1.grid(row=1, column=1)
label2 = Label(window, text='Aerofoil name', width=label_width, anchor=W)
label2.grid(row=2, column=1)
combo = Combobox(window, width=text_width - 3)
combo['values'] = ('NACA0006', 'NACA0008', 'NACA0009', 'NACA0010', 'NACA0012', 'NACA0015', 'NACA0018', 'NACA0021',
                   'NACA0024', 'NACA1408', 'NACA1410', 'NACA1412', 'NACA2408', 'NACA2410', 'NACA2411', 'NACA2412',
                   'NACA2414', 'NACA2415', 'NACA2418', 'NACA2421', 'NACA2424', 'NACA4412', 'NACA4415', 'NACA4418',
                   'NACA4421', 'NACA4424', 'NACA6409', 'NACA6412')
combo.current(15)
combo.grid(row=2, column=2)
# text2 = Entry(window, width=text_width)
# text2.grid(row=2, column=2)

label3 = Label(window, text='Flow parameters', font=('Arial Black', 10), width=label_width, anchor=W)
label3.grid(row=3, column=1)

label4 = Label(window, text='Reynolds number', width=label_width, anchor=W)
label4.grid(row=4, column=1)
text4 = Entry(window, width=text_width)
text4.grid(row=4, column=2)
text4.insert(END, '1e6')

label5 = Label(window, text='Density', width=label_width, anchor=W)
label5.grid(row=5, column=1)
text5 = Entry(window, width=text_width)
text5.grid(row=5, column=2)
text5.insert(END, '1.225')

label6 = Label(window, text='Viscosity', width=label_width, anchor=W)
label6.grid(row=6, column=1)
text6 = Entry(window, width=text_width)
text6.grid(row=6, column=2)
text6.insert(END, '1.789e-5')

label7 = Label(window, text='Angle of attack - degrees', width=label_width, anchor=W)
label7.grid(row=7, column=1)
text7 = Entry(window, width=text_width)
text7.grid(row=7, column=2)
text7.insert(END, '0.0')

label8 = Label(window, text='Plunging parameters', font=('Arial Black', 10), width=label_width, anchor=W)
label8.grid(row=8, column=1)

label9 = Label(window, text='Plunging amplitude - m', width=label_width, anchor=W)
label9.grid(row=9, column=1)
text9 = Entry(window, width=text_width)
text9.grid(row=9, column=2)
text9.insert(END, '1.0')

label10 = Label(window, text='Plunging frequency - Hz', width=label_width, anchor=W)
label10.grid(row=10, column=1)
text10 = Entry(window, width=text_width)
text10.grid(row=10, column=2)
text10.insert(END, '5.0')

label11 = Label(window, text='Time parameters', font=('Arial Black', 10), width=label_width, anchor=W)
label11.grid(row=11, column=1)

label12 = Label(window, text='Time step - s', width=label_width, anchor=W)
label12.grid(row=12, column=1)
text12 = Entry(window, width=text_width)
text12.grid(row=12, column=2)
text12.insert(END, '0.005')

label13 = Label(window, text='Iteration', width=label_width, anchor=W)
label13.grid(row=13, column=1)
text13 = Entry(window, width=text_width)
text13.grid(row=13, column=2)
text13.insert(END, '1000')

label14 = Label(window, text='New vortex placement', font=('Arial Black', 10), width=label_width, anchor=W)
label14.grid(row=14, column=1)

label15 = Label(window, text='Distance - %chord', width=label_width, anchor=W)
label15.grid(row=15, column=1)
text15 = Entry(window, width=text_width)
text15.grid(row=15, column=2)
text15.insert(END, '0.001')

label16 = Label(window, text='Angle - degrees', width=label_width, anchor=W)
label16.grid(row=16, column=1)
text16 = Entry(window, width=text_width)
text16.grid(row=16, column=2)
text16.insert(END, '0.0')

# column 2
label01 = Label(window, text='         ')
label01.grid(row=1, column=3)

# column 3
label17 = Label(window, text='Data files', font=('Arial Black', 10), width=10, anchor=W)
label17.grid(row=1, column=4, sticky=W)

chk_state1 = BooleanVar()
chk_state1.set(True)
checkbtn1 = Checkbutton(window, text='Vortex position file', var=chk_state1)
checkbtn1.grid(row=2, column=4, sticky=W)
chk_state2 = BooleanVar()
chk_state2.set(True)
checkbtn2 = Checkbutton(window, text='Force calculation file', var=chk_state2)
checkbtn2.grid(row=3, column=4, sticky=W)
chk_state3 = BooleanVar()
chk_state3.set(False)
# checkbtn3 = Checkbutton(window, text='Excel file', var=chk_state3)
# checkbtn3.grid(row=4, column=4, sticky=W)

label18 = Label(window, text='Graph plot', font=('Arial Black', 10), width=10, anchor=W)
label18.grid(row=5, column=4, sticky=W)

chk_state4 = BooleanVar()
chk_state4.set(True)
checkbtn4 = Checkbutton(window, text='Vortex movement', var=chk_state4)
checkbtn4.grid(row=6, column=4, sticky=W)
chk_state5 = BooleanVar()
chk_state5.set(False)
checkbtn5 = Checkbutton(window, text='Cl vs time', var=chk_state5)
checkbtn5.grid(row=7, column=4, sticky=W)
chk_state6 = BooleanVar()
chk_state6.set(False)
checkbtn6 = Checkbutton(window, text='Cd vs time', var=chk_state6)
checkbtn6.grid(row=8, column=4, sticky=W)
chk_state7 = BooleanVar()
chk_state7.set(False)
checkbtn7 = Checkbutton(window, text='Plunging distance vs time', var=chk_state7)
checkbtn7.grid(row=9, column=4, sticky=W)

button1 = Button(window, text='Calculate', command=button_command)
button1.grid(row=17, column=4, sticky=E)
button2 = Button(window, text='Stop', command=exit_command)
button2.grid(row=18, column=4, sticky=E)

window.geometry('530x430')  # windwo size customization
window.resizable(0, 0)
window.mainloop()
