import numpy as np
from matplotlib import pyplot as plt

A = np.array([[1.6, -0.5,0,0,0], [1.1,-0.6,-0.5,0,0],[0,1.1,-0.6,-0.5,0],[0,0,1.1,-0.6,-0.5],[0,0,0,1.6,-0.6]])
B = np.array([1.1,0,0,0,0])
X = np.linalg.inv(A).dot(B)
print(X)
'''
# map into z plane
radius = 0.27546
center_circle = 0.48333 + 0.01141j
num_div = 100
phi = 2 * np.pi / num_div  # divisable number
angle = np.array([n * phi for n in range(num_div)])
circle_point_u = np.exp(1j * angle)
circle_point_v = circle_point_u * radius + center_circle

Gkn = np.array([2.2191e-01 - 1.4299e-03j, 1.5057e-02 - 8.7489e-03j, 1.4982e-03 - 1.1215e-04j,
                9.0373e-04 - 2.4252e-04j, 4.3968e-04 - 2.9267e-04j, 3.0882e-04 - 1.4200e-05j,
                1.4553e-04 - 2.4458e-05j, 1.3449e-04 - 5.7484e-05j, 6.2759e-05 - 4.8699e-05j,
                6.4421e-05 + 3.4589e-06j, 3.8072e-05 - 1.8774e-05j, 3.7293e-05 - 2.0058e-05j,
                2.1524e-05 - 1.3357e-05j, 2.2883e-05 + 3.4335e-07j, 1.5773e-05 - 1.1712e-05j,
                1.4867e-05 - 8.2966e-06j, 1.0312e-05 - 4.6961e-06j, 1.0592e-05 - 1.2576e-06j,
                8.1001e-06 - 7.2054e-06j, 7.3614e-06 - 3.5700e-06j, 6.0316e-06 - 2.1841e-06j,
                5.8593e-06 - 1.7983e-06j, 4.7887e-06 - 4.2143e-06j, 4.2513e-06 - 1.6065e-06j,
                3.8882e-06 - 1.3190e-06j, 3.7131e-06 - 1.7341e-06j, 3.0863e-06 - 2.4206e-06j,
                2.9155e-06 - 8.1604e-07j, 2.6536e-06 - 1.0478e-06j, 2.6585e-06 - 1.3842e-06j,
                2.1356e-06 - 1.3556e-06j, 2.1907e-06 - 5.5524e-07j, 1.9991e-06 - 9.2462e-07j,
                2.0653e-06 - 1.0256e-06j, 1.5499e-06 - 8.6178e-07j, 1.8207e-06 - 4.3190e-07j,
                1.4895e-06 - 8.5323e-07j, 1.6596e-06 - 7.0904e-07j, 1.1445e-06 - 4.9074e-07j,
                1.5170e-06 - 5.0824e-07j, 1.2971e-06 - 6.2508e-07j, 1.2989e-06 - 5.4237e-07j,
                9.8841e-07 - 4.3777e-07j, 1.4495e-06 - 6.7442e-07j, 1.0236e-06 - 2.3024e-07j,
                1.0254e-06 - 7.7694e-07j, 1.2409e-06 - 1.9772e-07j, 9.0645e-07 - 5.8776e-07j,
                1.0258e-06 - 3.9214e-07j, 1.0181e-06 - 4.7850e-07j])

# circle_point_u = np.array([1, 2, 3, 4])
# Gkn = np.array([10, 20, 30])
# radius = 1
# center_circle = 10000

n = np.arange(1, len(Gkn) + 1)
circle_point_u_1 = np.tile(circle_point_u, (len(Gkn), 1))
Gkn_1 = np.tile(Gkn, (len(circle_point_u), 1)).transpose()
n_1 = np.tile(n, (len(circle_point_u), 1)).transpose()

airfoil_point = circle_point_u * radius + center_circle + sum(Gkn_1 / (circle_point_u_1 ** n_1))
print(airfoil_point)

plt.plot(circle_point_u.real, circle_point_u.imag)
plt.gca().set_aspect('equal', adjustable='box')
plt.scatter(circle_point_u.real, circle_point_u.imag, s=2)

plt.plot(circle_point_v.real, circle_point_v.imag)
plt.gca().set_aspect('equal', adjustable='box')
plt.scatter(circle_point_v.real, circle_point_v.imag, s=2)

plt.scatter(airfoil_point.real, airfoil_point.imag, s=2)
plt.show()
'''
'''
# circle point
radius = 0.27546
center_circle = 0.48333+0.01141j
num_div = 50
phi = 2 * np.pi / num_div  # divisable number
angle = np.array([n*phi for n in range(num_div)])
circle_point_u = center_circle + radius * np.exp(1j * angle)
circle_point_v = circle_point_u * radius + center_circle

plt.plot(circle_point_u.real, circle_point_u.imag)
plt.gca().set_aspect('equal', adjustable='box')
plt.scatter(circle_point_u[0].real, circle_point_u[0].imag)

plt.plot(circle_point_v.real, circle_point_v.imag)
plt.gca().set_aspect('equal', adjustable='box')
plt.scatter(circle_point_v.real, circle_point_v.imag)
plt.show()
'''
