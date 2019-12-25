from matplotlib import pyplot as plt
import numpy as np


def plot_graph_all(axs, heading_list, x_axis_title, y_axis_title, x_data, y_data, pause_time, islast, steady_state_cl,
                   x_limit_high, x_limit_low, path_to_save, plunging_on):
    # iteration vs time
    axs[0, 0].set_xlabel(x_axis_title[0])
    axs[0, 0].set_ylabel(y_axis_title[0])
    axs[0, 0].set_title(heading_list[0])
    axs[0, 0].plot(x_data[0], y_data[0], color='B', label='plot 1', linewidth=1)
    axs[0, 0].grid(True)
    # cl/cd
    axs[0, 1].set_xlabel(x_axis_title[1])
    axs[0, 1].set_ylabel(y_axis_title[1])
    axs[0, 1].set_title('Cl &Cd')
    axs[0, 1].set_xlim([x_limit_low, x_limit_high])
    axs[0, 1].plot(x_data[1], y_data[1], color='R', label='plot 1', linewidth=1)
    axs[0, 1].plot(x_data[2], y_data[2], color='B', label='plot 1', linewidth=1)
    if not plunging_on:
        axs[0, 1].plot(x_data[1], steady_state_cl, color='y', label='plot 1', linewidth=1)
    axs[0, 1].grid(True)
    # cl
    axs[1, 0].set_xlabel(x_axis_title[1])
    axs[1, 0].set_ylabel(y_axis_title[1])
    axs[1, 0].set_title(heading_list[1])
    axs[1, 0].set_xlim([x_limit_low, x_limit_high])
    axs[1, 0].plot(x_data[1], y_data[1], color='R', label='plot 1', linewidth=1)
    if not plunging_on:
        axs[1, 0].plot(x_data[1], steady_state_cl, color='y', label='plot 1', linewidth=1)
    axs[1, 0].grid(True)
    # cd
    axs[1, 1].set_xlabel(x_axis_title[2])
    axs[1, 1].set_ylabel(y_axis_title[2])
    axs[1, 1].set_title(heading_list[2])
    axs[1, 1].set_xlim([x_limit_low, x_limit_high])
    axs[1, 1].plot(x_data[2], y_data[2], color='B', label='Cl', linewidth=1)
    axs[1, 1].grid(True)

    plt.tight_layout()
    plt.pause(pause_time)

    if islast:
        plt.savefig(path_to_save + '/Figure.png', dpi=1280)
        plt.show()


def cl_grapgh_plot(x1, y1, pause_time, islast, x_limit_low, x_limit_high, steady, path_to_save, plunging_on):
    plt.figure(1)
    plt.clf()
    plt.xlabel('Time (s)')
    plt.ylabel('Cl')
    plt.title('Cl vs time')
    plt.xlim(x_limit_low, x_limit_high)
    plt.plot(x1, y1, color='r', label='transient cl', linewidth=1)
    if not plunging_on:
        plt.plot(x1, steady, color='b', label='steady state cl', linewidth=1)
    plt.grid(True)
    plt.tight_layout()
    plt.pause(pause_time)
    plt.legend()
    if islast:
        plt.savefig(path_to_save + '/Cl vs time.png')
        plt.show()


def cd_grapgh_plot(x1, y1, pause_time, islast, x_limit_low, x_limit_high, path_to_save):
    plt.figure(2)
    plt.clf()
    plt.xlabel('Time (s)')
    plt.ylabel('Cd')
    plt.title('Cd vs time')
    plt.xlim(x_limit_low, x_limit_high)
    plt.plot(x1, y1, color='r', linewidth=1)
    plt.grid(True, which='both')
    plt.tight_layout()
    plt.pause(pause_time)

    if islast:
        plt.savefig(path_to_save + '/Cd vs time.png')
        plt.show()


def plunging_distance_plot(x1, y1, pause_time, islast, x_limit_low, x_limit_high, path_to_save):
    plt.figure(3)
    plt.clf()
    plt.xlabel('Time (s)')
    plt.ylabel('distance')
    plt.title('Plunging wing movement')
    plt.xlim(x_limit_low, x_limit_high)
    plt.plot(x1, y1, color='R', linewidth=1)
    plt.grid(True, which='both')
    plt.tight_layout()
    plt.pause(pause_time)

    if islast:
        plt.savefig(path_to_save + '/plunging movement.png')
        plt.show()


def plot_airfoil(airfoil, vortex, strength, aoa, pause_time, islast):  # free stream initial aoa should there
    aoa = np.exp(-1j * aoa)
    vortex_pos = np.array([])
    vortex_neg = np.array([])
    for index in range(len(vortex)):
        if strength[index] < 0:
            vortex_neg = np.append(vortex_neg, [vortex[index]])
        else:
            vortex_pos = np.append(vortex_pos, [vortex[index]])
    vortex_pos = vortex_pos * aoa
    vortex_neg = vortex_neg * aoa
    airfoil = airfoil * aoa

    plt.figure(4)
    plt.clf()
    plt.axis('off')
    plt.title('vortex movement')
    plt.grid(False)
    plt.axis('equal')
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(vortex_pos.real, vortex_pos.imag, s=2, color='b')
    plt.scatter(vortex_neg.real, vortex_neg.imag, s=2, color='r')
    plt.plot(airfoil.real, airfoil.imag)

    plt.tight_layout()
    plt.pause(pause_time)
    if islast:
        # plt.savefig(path_to_save + '/plunging movement.png')
        plt.show()


'''
plot parameters
plt.plot(dev_x, dev_y, label='plot 1', linewidth=1)
'''

'''
scatter parameters
plt.scatter(dev_x, dev_y2, s=50, c=dev_y2, cmap='Reds', label='name')
cbar = plt.colorbar()
cbar.set_label('satisfaction')
'''
