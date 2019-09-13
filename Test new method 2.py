import function_file as ff
import numpy as np
import time


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
free_aoa = np.deg2rad(free_aoa)
# ------ plunging parameters
pl_amplitude = 0
pl_frequency = 0
# ------ pitching parameters
pi_amplitude = 0
pi_frequency = 0
# ------ new vortex
distance = 0.005
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
time_step = 0.01
current_time = 0.00
iteration = 1

# ------ derivatives

# ----- write in a file
# ff.make_file(airfoil, N, radius, center_circle, trailing_edge_z, Gkn, z_plane, v_plane, u_plane,
#              free_velocity, free_aoa, pl_amplitude, pl_frequency, pi_amplitude, pi_frequency,
#              time_step, current_time, iteration, distance, angle)

# ------ iteration code
for iterate in range(iteration):
    a = np.array(range(len(Gkn))) + 1

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
    te_vortex_vel = np.conj(vel_conj)
    te_vortex_vel = sum(np.transpose(te_vortex_vel))
    
print('total time ', time.time() - start)
