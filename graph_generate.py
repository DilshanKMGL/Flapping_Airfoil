import matplotlib.pyplot as plt
import numpy as np
import function_file as ff

file1 = open('result_file.txt', 'r')
line = file1.readlines()

center = complex(line[5][1:-2])
iteration = int(line[33][:-1])
airfoil_z = line[11][1:-2].split(',')
airfoil_z = [complex(coord) for coord in airfoil_z]
airfoil_z_real = [index.real for index in airfoil_z]
airfoil_z_imag = [index.imag for index in airfoil_z]
heading = str(line[1][:-2])
te_vortex_strength_line = 44
te_vortex_z_line = 46

for i in range(iteration):
    print('Iteration :'+str(i))
    te_vortex_strength = line[te_vortex_strength_line][1:-2]
    te_vortex_strength = [float(index) for index in te_vortex_strength.split(',')]

    te_vortex_z = line[te_vortex_z_line][1:-2].replace('*I', 'j').replace(' ', '')
    te_vortex_z = [complex(index) for index in te_vortex_z.split(',')]

    te_vortex_real_part_positive = []
    te_vortex_imag_part_positive = []
    te_vortex_real_part_negative = []
    te_vortex_imag_part_negative = []

    for index in range(len(te_vortex_strength)):
        if te_vortex_strength[index] > 0:
            te_vortex_real_part_positive.append(te_vortex_z[index].real)
            te_vortex_imag_part_positive.append(te_vortex_z[index].imag)
        else:
            te_vortex_real_part_negative.append(te_vortex_z[index].real)
            te_vortex_imag_part_negative.append(te_vortex_z[index].imag)

    te_vortex_strength_line += 5
    te_vortex_z_line += 5

    plt.axis('off')
    #plt.xlim(-1, 10)
    #plt.ylim(-5, 5)

    plt.grid(False)
    plt.plot(airfoil_z_real, airfoil_z_imag, color='b')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(te_vortex_real_part_positive, te_vortex_imag_part_positive, s=2, color='g')
    plt.scatter(te_vortex_real_part_negative, te_vortex_imag_part_negative, s=2, color='r')

    temp_heading = heading + ' ' + str(i)
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.savefig(temp_heading)
    plt.close()
'''
time_step = float(line[6][21:])
pl_amp = float(line[9][11:])
pl_fre = float(line[10][11:])
iterate_number = int(line[7][21:])
strength_list = line[len(line) - 7][1:len(line[len(line) - 7]) - 2].split(',')
strength = [float(index) for index in strength_list]

# ---------------------- airfoil --------------------------
r = np.sqrt(center.imag ** 2 + (a - center.real) ** 2)
theta = np.linspace(0, 2.0 * np.pi, 20)  # 200000
cir = center + r * np.exp(1j * theta)
cir = np.array(cir)
jou = cir + a ** 2 / cir

# ---------------------- plot -----------------------------
vortex_iteration = 18
vortex_line = 21

for i in range(iterate_number):
    iteration = int(line[vortex_iteration][10:])
    print("iteration ", iteration)
    func = ff.plunging(pl_amp, pl_fre, time_step * (iteration - 1))
    jou_new = jou + complex(str(func[0]).replace('*', '').replace('I', 'j'))

    vor_list = line[vortex_line][1:len(line[vortex_line]) - 2].split(',')
    vortex = [complex(index.replace(' ', '').replace('I', 'j').replace('*', '')) for index in vor_list]

    real_part_positive = []
    imag_part_positive = []
    real_part_negative = []
    imag_part_negative = []
    for index in range(len(vortex)):
        if strength[index] > 0:
            real_part_positive.append(vortex[index].real)
            imag_part_positive.append(vortex[index].imag)
        else:
            real_part_negative.append(vortex[index].real)
            imag_part_negative.append(vortex[index].imag)
    plt.axis('off')
    plt.xlim(-5, 60)
    plt.ylim(-20, 20)

    plt.grid(False)
    plt.plot(jou_new.real, jou_new.imag, color='b')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(real_part_positive, imag_part_positive, s=2, color='g')
    plt.scatter(real_part_negative, imag_part_negative, s=2, color='r')

    heading = 'Plunging pattern : ' + str(iteration)
    plt.savefig(heading)
    plt.close()

    vortex_iteration += 5
    vortex_line += 5
'''
