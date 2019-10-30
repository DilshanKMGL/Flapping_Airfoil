import numpy as np


def diag_remove(p):
    strided = np.lib.stride_tricks.as_strided

    p = np.tile(p, (len(p), 1))
    m = p.shape[0]
    s0, s1 = p.strides
    p = strided(p.ravel()[1:], shape=(m - 1, m), strides=(s0 + s1, s1)).reshape(m, -1).transpose()
    return p


def calculate_derivative(point):  # calculate derivative of a single point
    dzdu = radius - sum(power*Gkn/point**(power+1))
    d2zdu2 = sum(power*(power+1)/point**(power+2))
    dudz = 1/dzdu
    d2udz2 = -d2zdu2/dzdu**3
    return dudz, d2udz2


def calculate_d4(point, point_aff, v_strength):  # calculate vortex induced velocity
    p1 = (point * radius + center_circle) + sum(Gkn / point ** power) - point_aff
    p2 = radius / point + sum(Gkn * (point * radius / (radius - point * center_circle)) ** power) - point_aff
    p3 = radius / point + sum(Gkn * (point * radius / (radius - point * center_circle)) ** power) - point

    p11 = radius - sum(power * Gkn / point ** (power + 1.0))
    p22 = -radius / point ** 2.0 + \
         sum(Gkn * power * point ** (power - 1) * (radius / (radius - point * center_circle)) ** (power + 1))
    p33 = p22

    d4 = v_strength * 1j / (np.pi * 2) * (p22 / p2 + p33 / p3 - p11 / p1)
    return d4


velocity = 20
aoa = 0
te_vortex_strength = np.array([2, 4])
te_vortex_u = np.array([30, 40])
Gkn = np.array([0.22191 - 0.0014299j, 0.015057 - 0.0087489j, 0.0014982 - 0.00011215j, 0.00090373 - 0.00024252j,
                0.00043968 - 0.00029267j, 0.00030882 - 1.42e-05j, 0.00014553 - 2.4458e-05j, 0.00013449 - 5.7484e-05j,
                6.2759e-05 - 4.8699e-05j, 6.4421e-05 + 3.4589e-06j, 3.8072e-05 - 1.8774e-05j, 3.7293e-05 - 2.0058e-05j,
                2.1524e-05 - 1.3357e-05j, 2.2883e-05 + 3.4335e-07j, 1.5773e-05 - 1.1712e-05j, 1.4867e-05 - 8.2966e-06j,
                1.0312e-05 - 4.6961e-06j, 1.0592e-05 - 1.2576e-06j, 8.1001e-06 - 7.2054e-06j, 7.3614e-06 - 3.57e-06j,
                6.0316e-06 - 2.1841e-06j, 5.8593e-06 - 1.7983e-06j, 4.7887e-06 - 4.2143e-06j, 4.2513e-06 - 1.6065e-06j,
                3.8882e-06 - 1.319e-06j, 3.7131e-06 - 1.7341e-06j, 3.0863e-06 - 2.4206e-06j, 2.9155e-06 - 8.1604e-07j,
                2.6536e-06 - 1.0478e-06j, 2.6585e-06 - 1.3842e-06j, 2.1356e-06 - 1.3556e-06j, 2.1907e-06 - 5.5524e-07j,
                1.9991e-06 - 9.2462e-07j, 2.0653e-06 - 1.0256e-06j, 1.5499e-06 - 8.6178e-07j, 1.8207e-06 - 4.319e-07j,
                1.4895e-06 - 8.5323e-07j, 1.6596e-06 - 7.0904e-07j, 1.1445e-06 - 4.9074e-07j, 1.517e-06 - 5.0824e-07j,
                1.2971e-06 - 6.2508e-07j, 1.2989e-06 - 5.4237e-07j, 9.8841e-07 - 4.3777e-07j, 1.4495e-06 - 6.7442e-07j,
                1.0236e-06 - 2.3024e-07j, 1.0254e-06 - 7.7694e-07j, 1.2409e-06 - 1.9772e-07j, 9.0645e-07 - 5.8776e-07j,
                1.0258e-06 - 3.9214e-07j, 1.0181e-06 - 4.785e-07j])
power = np.arange(1, len(Gkn) + 1, dtype='float')
radius = 0.27546
center_circle = 0.48333 + 0.01141j

vortex_strength_1 = 50
vortex_strength_2 = 30
vortex_center_1 = 3050 + 5j
vortex_center_2 = 3050 - 5j

# volocities in z plane - calulate in u plane
point = vortex_center_1
point_aff = vortex_center_2
d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / point ** 2)
dudz, d2udz2 = calculate_derivative(point)

p1 = (point * radius + center_circle) + sum(Gkn / point ** power) - point_aff
p11 = radius - sum(power * Gkn / point ** (power + 1.0))
p2 = (radius / point) + center_circle + sum(Gkn * point ** power) - np.conjugate(point_aff)
p22 = -radius / point ** 2.0 + sum(Gkn * power * point ** (power - 1))
p3 = (radius / point) + center_circle + sum(Gkn * point ** power) - np.conjugate(point)
p33 = p22

d4 = 1j / (2 * np.pi) * vortex_strength_2 * (p22/p2 - p11/p1) + 1j / (2 * np.pi) * vortex_strength_1 * (p33/p3)
p = d4
vel_1 = np.conjugate(p*dudz - 1j*vortex_strength_1 / (4*np.pi) * d2udz2/dudz)
print('veolcity vortex1 - uplane ', vel_1)

# velocities in z plane - calculate in z plane
vortex_center_1_v = vortex_center_1 * radius + center_circle
vortex_center_2_v = vortex_center_2 * radius + center_circle
vortex_center_1_z = vortex_center_1_v + sum(Gkn * (radius / (vortex_center_1_v - center_circle)) ** power)
vortex_center_2_z = vortex_center_2_v + sum(Gkn * (radius / (vortex_center_2_v - center_circle)) ** power)

# only 2 vortices
vel_1_z = - 1j*vortex_strength_2 / (2 * np.pi) / (vortex_center_1_z - vortex_center_2_z)
print('veolcity vortex1 - zplane ', vel_1_z)

# one vortex + free stream
# vel_1_z = velocity
# print('veolcity vortex1 - zplane ', vel_1_z)