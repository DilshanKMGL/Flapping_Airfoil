import numpy as np
from matplotlib import pyplot as plt

center = -0.05 + 1j * 0.05
a = 0.5
r = np.sqrt(center.imag ** 2 + (a - center.real) ** 2)
theta = np.linspace(0, 2.0 * np.pi, 100)  # 200000
cir = center + r * np.exp(1j * theta)
cir = np.array(cir)
jou = cir + a ** 2 / cir

x = jou.real
y = jou.imag

file1 = open('joukowski_airfoil.txt', 'w')
file1.write('[')
for index in range(len(jou)):
    file1.write(str(jou.real[index])+'+'+str(jou.imag[index])+'i\n')
file1.write(']')
file1.close()


plt.grid(False)
plt.plot(jou.real, jou.imag, color='b')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
