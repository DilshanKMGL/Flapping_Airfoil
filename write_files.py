import numpy as np
# import xlsxwriter as xl


def make_file(airfoil, re_num, strauhl_number, density, viscosity, free_velocity, free_aoa, pl_amplitude, pl_frequency,
              pi_amplitude, pi_frequency, time_step, current_time, iteration, distance, angle, heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, 'w')

    file1.write('airfoil\n' + str(airfoil) + '\n')
    file1.write('re_num\n' + str(re_num) + '\n')
    file1.write('strauhl number\n' + str(strauhl_number) + '\n')
    file1.write('density\n' + str(density) + '\n')
    file1.write('viscosity\n' + str(viscosity) + '\n')
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
    file1.write('te_vortex_z' + '\n')
    file1.write('First value is the iteration number' + '\n')
    file1.close()


def make_force_file(airfoil, re_num, strauhl_number, density, viscosity, free_velocity, free_aoa, pl_amplitude,
                    pl_frequency, pi_amplitude, pi_frequency, time_step, current_time, iteration, distance, angle,
                    heading):
    # heading = 'Transient_solution_results/' + 'result_file_' + airfoil + '.txt'
    file1 = open(heading, 'w')

    file1.write('airfoil\n' + str(airfoil) + '\n')
    file1.write('re_num\n' + str(re_num) + '\n')
    file1.write('strauhl number\n' + str(strauhl_number) + '\n')
    file1.write('density\n' + str(density) + '\n')
    file1.write('viscosity\n' + str(viscosity) + '\n')
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


def create_excel_file(path_dir):
    force_file = open(path_dir + '/force_file.txt')
    result_file = open(path_dir + '/result_file.txt')
    force_line = force_file.readlines()
    result_line = result_file.readlines()
    airfoil_name = force_line[1][:-1]
    # workbook = xl.Workbook(path_dir + '/Data.xlsx')
    # worksheet = workbook.add_worksheet(airfoil_name)
    # workbook.close()
