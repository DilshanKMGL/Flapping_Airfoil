import function_file as ff
import sympy as sp
import numpy as np
import time
import cmath


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
# ------ new vortex
distance = 0.005
angle = 0
# ------ data store
circulation_list = np.array([])
te_vortex_strength = np.array([])
te_vortex_z = np.array([])
te_vortex_u = np.array([])
iterate_time_step = np.array([])
# ------ time step
time_step = 0.01
current_time = 0.00
iteration = 1

# ------ create functions
v, u, vor = sp.symbols('v, u, vor', real=False)
n = np.arange(1, len(Gkn) + 1)
z_fun = v + sum(Gkn * (radius / (v - center_circle)) ** n)
v_fun = z_fun.subs(v, u * radius + center_circle)
# ------ derivatives

# ----- write in a file
# ff.make_file(airfoil, N, radius, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane,
#              free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
#              time_step, current_time, iteration, distance, angle)

# ------ iteration code
for iterate in range(iteration):
    a = np.array(range(len(Gkn)))+1

    print('Iteration - ' + str(iterate + 1))

    # ------ calculate trailing edge position
    search_point = center_circle + radius
    trailing_edge_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, trailing_edge_z))
    trailing_edge_u = complex((trailing_edge_v - center_circle) / radius)
    # --- calculate velocity
    velocity = free_velocity
    aoa = free_aoa

    # ------ calculate new vortex position
    new_vortex_position_z = trailing_edge_z + distance * sp.exp(-1j * sp.rad(angle)).evalf()
    new_vortex_position_z = complex(new_vortex_position_z)
    search_point = center_circle + radius
    new_vortex_position_v = complex(newton(search_point, 1e-8, 50, Gkn, radius, center_circle, new_vortex_position_z))
    new_vortex_position_u = complex((new_vortex_position_v - center_circle) / radius)
    te_vortex_z = np.append(te_vortex_z, [new_vortex_position_z])
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
    circulation = - (s * d2 + d3 + d4) / (d1 + d2)

    print(circulation.evalf())
    print('Iteratation ' + str(iterate + 1) + ' complete after ')

print('total time ', round(time.time() - start, 2))
