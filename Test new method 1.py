import function_file as ff
import sympy as sp
import numpy as np
import time
import cmath


def newton(f, Df, x0, epsilon, max_iter):
    """Approximate solution of f(x)=0 by Newton's method.

    Parameters
    ----------
    f : function
        Function for which we are searching for a solution f(x)=0.
    Df : function
        Derivative of f(x).
    x0 : number
        Initial guess for a solution f(x)=0.
    epsilon : number
        Stopping criteria is abs(f(x)) < epsilon.
    max_iter : integer
        Maximum number of iterations of Newton's method.

    Returns
    -------
    xn : number
        Implement Newton's method: compute the linear approximation
        of f(x) at xn and find x intercept by the formula
            x = xn - f(xn)/Df(xn)
        Continue until abs(f(xn)) < epsilon and return xn.
        If Df(xn) == 0, return None. If the number of iterations
        exceeds max_iter, then return None.

    Examples
    --------
    #>>> f = lambda x: x**2 - x - 1
    #>>> Df = lambda x: 2*x - 1
    #>>> newton(f,Df,1,1e-8,10)
    Found solution after 5 iterations.
    1.618033988749989
    """

    v = sp.symbols('v', real=False)
    xn = x0
    for n in range(0, max_iter):
        fxn = f.subs([(v, xn)]).evalf()  # f(xn)

        if abs(fxn) < epsilon:
            # print('Found solution after', n, 'iterations.')
            return xn
        Dfxn = Df.subs([(v, xn)]).evalf()  # Df(xn)
        if Dfxn == 0:
            # print('Zero derivative. No solution found.')
            return None
        xn = (xn - fxn / Dfxn).evalf()
    print('Exceeded maximum iterations. No solution found.')
    return None


start = time.time()
iterate_time = start
# ------ airfoil data
airfoil = 'NACA2412'
N, radius, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane = ff.read_data(airfoil)
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
iteration = 25
# ------ new vortex
distance = 0.005
angle = 0
# ------ data store
circulation = np.array([])
te_vortex_strength = np.array([])
te_vortex_z = np.array([])
te_vortex_u = np.array([])
iterate_time_step = np.array([])

# ------ create functions
v, u, vor = sp.symbols('v, u, vor', real=False)
n = np.arange(1, len(Gkn) + 1)
z_fun = v + sum(Gkn * (radius / (v - center_circle)) ** n)
v_fun = z_fun.subs(v, u * radius + center_circle)
# ------ derivatives
first_derivative = 1 / sp.diff(v_fun, u)
second_derivative = sp.diff(1 / sp.diff(v_fun, u), u) / sp.diff(v_fun, u)
# ----- write in a file
ff.make_file(airfoil, N, radius, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane,
             free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
             time_step, current_time, iteration, distance, angle)

# ------ iteration code
for iterate in range(iteration):
    print('Iteration - ' + str(iterate + 1))
    # ------ calculate trailing edge position
    search_point = center_circle + radius
    trailing_edge_v = newton(z_fun - trailing_edge_z, sp.diff(z_fun, v), search_point, 1e-8, 50)
    trailing_edge_u = (trailing_edge_v - center_circle) / radius

    # ------ calculate new vortex position
    new_vortex_position_z = trailing_edge_z + distance * sp.exp(-1j * sp.rad(angle)).evalf()
    search_point = center_circle + radius
    new_vortex_position_v = newton(z_fun - new_vortex_position_z, sp.diff(z_fun, v), search_point, 1e-8, 50)
    new_vortex_position_u = (new_vortex_position_v - center_circle) / radius
    te_vortex_z = np.append(te_vortex_z, [new_vortex_position_z])
    te_vortex_u = np.append(te_vortex_u, [new_vortex_position_u])

    # ------ create function to calculate circulation
    # - circulation
    init_circulation = -1j * vor * sp.log(u) / (2 * sp.pi)
    # - new vortex
    te_vortex_sum = sum(te_vortex_strength)
    func1 = -1j * (- te_vortex_sum - vor) * sp.log(u - te_vortex_u[-1]) / (2 * sp.pi)
    func2 = 1j * (- te_vortex_sum - vor) * sp.log((1 / u) - te_vortex_u[-1]) / (2 * sp.pi)

    # - freestream
    func3 = free_velocity * u * pow(np.e, -1j * free_aoa)
    func4 = free_velocity * (1 / u) * pow(np.e, 1j * free_aoa)

    # - shed vortex
    temp_func = 1j * te_vortex_strength * np.array([sp.log(((1 / u) - te_vortex_u[index]) / (u - te_vortex_u[index]))
                                                    for index in range(len(te_vortex_strength))]) / (2 * sp.pi)
    func5 = sum(temp_func)
    # - solve equation
    vortex_function = init_circulation + func1 + func2 + func3 + func4 + func5
    func = sp.diff(vortex_function, u)
    func = func.subs([(u, trailing_edge_u)])
    circulation_new = sp.solve(func)
    # - add values to the list
    circulation_new = complex(circulation_new[0]).real
    circulation = np.append(circulation, [circulation_new])
    te_vortex_strength = np.append(te_vortex_strength, [-te_vortex_sum - circulation_new])
    # print(te_vortex_strength) ###################################################################333
    # - iteration time
    iterate_time_step = np.append(iterate_time_step, [round(time.time() - iterate_time, 2)])
    iterate_time = time.time()
    # - write iteration step time in file
    ff.write_array(circulation, te_vortex_strength, te_vortex_u, te_vortex_z, iterate_time_step, iterate)
    # - increment time
    current_time += time_step
    current_time = round(current_time, 2)

    # ------ move vortex
    # - newly shed vortex
    func_vel_1 = -1j * te_vortex_strength[-1] * sp.log(u - te_vortex_u[-1]) / (2 * sp.pi)
    func_vel_2 = 1j * te_vortex_strength[-1] * sp.log((1 / u) - te_vortex_u[-1]) / (2 * sp.pi)
    # - circulation
    func_vel_3 = -1j * circulation_new * sp.log(u) / (2 * sp.pi)

    # func3, func4 - freestream
    # func5 - sum(temp_func)
    vortex_function = func3 + func4 + func5 + func_vel_1 + func_vel_2 + func_vel_3
    temp_func = np.append(temp_func, [func_vel_1 + func_vel_2])
    vortex_function = vortex_function - temp_func

    velocity_function = np.array([sp.diff(index, u) for index in vortex_function])
    velocity_function = velocity_function * first_derivative - \
                        (te_vortex_strength * 1j / sp.pi * second_derivative / first_derivative)
    velocity_function = [velocity_function[index].subs(u, te_vortex_u[index]).evalf()
                          for index in range(len(velocity_function))]

    velocity_function = np.conj(velocity_function)
    past_te_vortex_u = te_vortex_u.copy()
    te_vortex_z = np.array(te_vortex_z + velocity_function * time_step)
    te_vortex_v = [newton(z_fun - te_vortex_z[index], sp.diff(z_fun, v),
                          past_te_vortex_u[index] * radius + center_circle, 1e-8, 50)
                       for index in range(len(te_vortex_z))]
    te_vortex_u = (np.array(te_vortex_v) - center_circle) / radius

    if iterate == iteration - 1:
        ff.final_position(te_vortex_z, te_vortex_u)

    print('Iteratation ' + str(iterate + 1) + ' complete after ', round(iterate_time_step[-1], 2), 's')

print('total time ', round(time.time() - start, 2))
