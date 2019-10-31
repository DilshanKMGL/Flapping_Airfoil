import numpy as np


def diag_remove(p):
    strided = np.lib.stride_tricks.as_strided

    p = np.tile(p, (len(p), 1))
    m = p.shape[0]
    s0, s1 = p.strides
    p = strided(p.ravel()[1:], shape=(m - 1, m), strides=(s0 + s1, s1)).reshape(m, -1).transpose()
    return p


def calculate_derivative(point):  # calculate derivative of a single point
    dzdu = radius - sum(power * Gkn / point ** (power + 1))
    d2zdu2 = sum(power * (power + 1) / point ** (power + 2))
    dudz = 1 / dzdu
    d2udz2 = -d2zdu2 / dzdu ** 3
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


Gkn = np.array([0.22191364662 - 0.00142990992j,
                0.01505722598 - 0.00874889365j,
                0.00149823803 - 0.00011215076j,
                0.00090373148 - 0.00024251526j,
                0.00043968377 - 0.00029267360j,
                0.00030882255 - 0.00001419958j,
                0.00014553022 - 0.00002445814j,
                0.00013448707 - 0.00005748407j,
                0.00006275852 - 0.00004869871j,
                0.00006442110 + 0.00000345894j,
                0.00003807247 - 0.00001877369j,
                0.00003729343 - 0.00002005799j,
                0.00002152385 - 0.00001335710j,
                0.00002288283 + 0.00000034335j,
                0.00001577258 - 0.00001171182j,
                0.00001486676 - 0.00000829661j,
                0.00001031177 - 0.00000469612j,
                0.00001059233 - 0.00000125761j,
                0.00000810013 - 0.00000720538j,
                0.00000736145 - 0.00000356997j,
                0.00000603161 - 0.00000218410j,
                0.00000585932 - 0.00000179831j,
                0.00000478871 - 0.00000421427j,
                0.00000425127 - 0.00000160646j,
                0.00000388823 - 0.00000131899j,
                0.00000371313 - 0.00000173414j,
                0.00000308634 - 0.00000242062j,
                0.00000291552 - 0.00000081604j,
                0.00000265362 - 0.00000104776j,
                0.00000265852 - 0.00000138416j,
                0.00000213562 - 0.00000135560j,
                0.00000219070 - 0.00000055524j,
                0.00000199909 - 0.00000092462j,
                0.00000206529 - 0.00000102563j,
                0.00000154990 - 0.00000086178j,
                0.00000182070 - 0.00000043190j,
                0.00000148947 - 0.00000085323j,
                0.00000165956 - 0.00000070904j,
                0.00000114451 - 0.00000049074j,
                0.00000151699 - 0.00000050824j,
                0.00000129712 - 0.00000062508j,
                0.00000129890 - 0.00000054237j,
                0.00000098841 - 0.00000043777j,
                0.00000144950 - 0.00000067442j,
                0.00000102357 - 0.00000023024j,
                0.00000102544 - 0.00000077694j,
                0.00000124092 - 0.00000019772j,
                0.00000090645 - 0.00000058776j,
                0.00000102579 - 0.00000039214j,
                0.00000101806 - 0.00000047850j])
power = np.arange(1, len(Gkn) + 1, dtype='float')
radius = 0.27546
center_circle = 0.483329 + 0.011410j
velocity = 0
aoa = 0

vortex_strength_1 = 10
vortex_strength_2 = 10
vortex_center_1 = 3000
vortex_center_2 = 3050

# free stream
d1 = velocity * radius * (np.exp(-1j * aoa) - np.exp(1j * aoa) / vortex_center_1 ** 2)

# in general method that was wrong before
p1 = 1 / (vortex_center_1 - vortex_center_2)
p2 = 1 / vortex_center_1 / (1 - vortex_center_1 * np.conj(vortex_center_2))
p = -1j * vortex_strength_2 / (2 * np.pi) * (p1 + p2)

p2 = 1 / vortex_center_1 / (1 - vortex_center_1 * np.conj(vortex_center_1))
p += -1j * vortex_strength_1 / (2 * np.pi) * p2
p += d1

# derivatives
dzdu = radius - sum(Gkn * power / vortex_center_1 ** (power + 1))
d2zdu2 = sum(Gkn * power * (power + 1) / vortex_center_1 ** (power + 2))
dudz = 1 / dzdu
d2udz2 = - d2zdu2 / dzdu ** 3

vel_1 = np.conjugate(p * dudz - 1j * vortex_strength_1 / (4 * np.pi) * d2udz2 / dudz)
print('veolcity vortex1 - uplane ', vel_1)

p1 = 1 / (vortex_center_2 - vortex_center_1)
p2 = 1 / vortex_center_2 / (1 - vortex_center_2 * np.conj(vortex_center_1))
p = -1j * vortex_strength_1 / (2 * np.pi) * (p1 + p2)

p2 = 1 / vortex_center_2 / (1 - vortex_center_2 * np.conj(vortex_center_2))
p += -1j * vortex_strength_2 / (2 * np.pi) * p2
p += d1

# derivatives
dzdu = radius - sum(Gkn * power / vortex_center_2 ** (power + 1))
d2zdu2 = sum(Gkn * power * (power + 1) / vortex_center_2 ** (power + 2))
dudz = 1 / dzdu
d2udz2 = - d2zdu2 / dzdu ** 3

vel_1 = np.conjugate(p * dudz - 1j * vortex_strength_2 / (4 * np.pi) * d2udz2 / dudz)
print('veolcity vortex2 - uplane ', vel_1)

# velocities in z plane - calculate in z plane
vortex_center_1_v = vortex_center_1 * radius + center_circle
vortex_center_2_v = vortex_center_2 * radius + center_circle
vortex_center_1_z = vortex_center_1_v + sum(Gkn * (radius / (vortex_center_1_v - center_circle)) ** power)
vortex_center_2_z = vortex_center_2_v + sum(Gkn * (radius / (vortex_center_2_v - center_circle)) ** power)

# only 2 vortices
vel_1_z = np.conj(velocity - 1j * vortex_strength_2 / (2 * np.pi) / (vortex_center_1_z - vortex_center_2_z))
vel_2_z = np.conj(velocity - 1j * vortex_strength_1 / (2 * np.pi) / (vortex_center_2_z - vortex_center_1_z))
print('veolcity vortex1 - zplane ', vel_1_z)
print('veolcity vortex2 - zplane ', vel_2_z)
