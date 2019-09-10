import function_file as ff
import sympy as sp
import numpy as np
import time
import cmath

start = time.time()
iterate_time = start
# ------ airfoil data
airfoil = 'NACA2412'
N, r, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane = ff.read_data(airfoil)
# ------ free stream velocity
free_velocity = 5
free_aoa = 0.0
free_aoa = sp.rad(free_aoa).evalf()
# ------ plunging parameters
pl_amplitude = 0
pl_frequency = 0
# ------ pitching parameters
pi_amplitude = 0
pi_frequency = 0
# ------ time step
time_step = 0.01
current_time = 0.00
iteration = 5
# ------ new vortex
distance = 0.005
angle = 0
# ------ data store
circulation = []
te_vortex_strength = []
te_vortex_z = []
te_vortex_u = []
iterate_time_step = []

# ------ derivatives
v, u = sp.symbols('v, u', real=False)
z_fun, v_fun = ff.get_v_function(Gkn, r, center_circle)
z_fun, v_fun = ff.get_v_function(Gkn, r, center_circle)
first_derivative = 1 / sp.diff(v_fun, u)
second_derivative = sp.diff(1 / sp.diff(v_fun, u), u) / sp.diff(v_fun, u)

# ----- write in a file
ff.make_file(airfoil, N, r, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane,
             free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
             time_step, current_time, iteration, distance, angle)
# ----- iteration code

for iterate in range(iteration):
    print('Iteration - ' + str(iterate + 1))

    fun_1 = ff.plunging(pl_amplitude, pl_frequency, current_time)
    func = cmath.polar(fun_1[1] + free_velocity)
    velocity = func[0]
    aoa = free_aoa + func[1]

    search_point = center_circle + r
    # te_value is a list [z value, u value]
    te_value = ff.get_trailing_edge(z_fun, trailing_edge_z, search_point, center_circle, r, fun_1[0])
    new_vortex_z, new_vortex_u = ff.new_vortex_position(z_fun, trailing_edge_z, search_point, angle, distance,
                                                        center_circle, r)

    te_vortex_u.append(new_vortex_u)
    te_vortex_z.append(new_vortex_z)

    # create function to make complex potential function
    vortex_function = ff.create_circulation(center_circle).subs(v, (u * r) / center_circle).evalf()
    te_vortex_sum = sum(te_vortex_strength)
    vortex_function += ff.create_vortex(te_vortex_u[-1], te_vortex_sum)
    vortex_function += ff.get_freestream(velocity, aoa)

    temp_func = [ff.get_vortex(te_vortex_u[vortex], te_vortex_strength[vortex])
                 for vortex in range(len(te_vortex_strength))]

    vortex_function += sum(temp_func)
    new_circulation = ff.calculate_circulation(vortex_function, te_value[1])
    circulation.append((sp.re(new_circulation)))
    te_vortex_strength.append((sp.re(-te_vortex_sum - new_circulation)))

    # ------ write in the file
    iterate_time_step.append(time.time() - iterate_time)

    iterate_time = time.time()
    ff.write_array(circulation, te_vortex_strength, te_vortex_u, te_vortex_z, iterate_time_step, iterate)

    current_time += time_step
    current_time = round(current_time, 2)

    # ------ move vortex
    vortex_function = sum(temp_func) + ff.get_circulation(circulation[-1]) + \
                      ff.get_vortex(te_vortex_u[-1], te_vortex_strength[-1]) + \
                      ff.get_freestream(velocity, aoa)
    vortex_function = sp.diff(vortex_function.evalf(), u)

    velocity_function = [(vortex_function -
                          sp.diff(ff.get_vortex(te_vortex_u[index], te_vortex_strength[index]), u)).evalf()
                         for index in range(len(te_vortex_strength))]

    velocity_function = np.multiply(velocity_function, first_derivative) - (
            np.multiply(te_vortex_strength, 1j) / sp.pi *
            second_derivative / first_derivative)

    velocity_function = [velocity_function[index].subs(u, te_vortex_u[index]).evalf()
                         for index in range(len(velocity_function))]
    velocity_function = np.conj(velocity_function)
    past_te_vortex_u = te_vortex_u.copy()
    te_vortex_z = te_vortex_z + velocity_function * time_step
    te_vortex_z = te_vortex_z.tolist()

    te_vortex_u = [ff.newton_u(v_fun - te_vortex_z[index], sp.diff(v_fun, u), past_te_vortex_u[index], 1e-8, 50)
                   for index in range(len(te_vortex_z))]

    if iterate == iteration - 1:
        ff.final_position(te_vortex_z, te_vortex_u)

    print('Iteratation ' + str(iterate + 1) + ' complete after ', round(iterate_time_step[-1], 2), 's')

print('total time ', round(time.time() - start, 2))
