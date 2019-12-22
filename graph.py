from matplotlib import pyplot as plt
import numpy as np


def plot_graph_all(axs, heading_list, x_axis_title, y_axis_title, x_data, y_data, pause_time, islast, steady_state_cl,
                   x_limit_high, x_limit_low):
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
    axs[0, 1].plot(x_data[1], steady_state_cl, color='y', label='plot 1', linewidth=1)
    axs[0, 1].grid(True)
    # cl
    axs[1, 0].set_xlabel(x_axis_title[1])
    axs[1, 0].set_ylabel(y_axis_title[1])
    axs[1, 0].set_title(heading_list[1])
    axs[1, 0].set_xlim([x_limit_low, x_limit_high])
    axs[1, 0].plot(x_data[1], y_data[1], color='R', label='plot 1', linewidth=1)
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
        plt.show()


def cl_grapgh_plot(x1, y1, pause_time, islast, x_limit_low, x_limit_high, steady):
    plt.figure(1)
    plt.clf()
    plt.xlabel('Time (s)')
    plt.ylabel('Cl')
    plt.title('Cl vs time')
    plt.xlim(x_limit_low, x_limit_high)
    plt.plot(x1, y1, color='R', linewidth=1)
    plt.plot(x1, steady, color='B', linewidth=1)
    plt.grid(True)
    plt.tight_layout()
    plt.pause(pause_time)
    if islast:
        plt.show()


def cd_grapgh_plot(x1, y1, pause_time, islast, x_limit_low, x_limit_high):
    plt.figure(2)
    plt.clf()
    plt.xlabel('Time (s)')
    plt.ylabel('Cd')
    plt.title('Cd vs time')
    plt.xlim(x_limit_low, x_limit_high)
    plt.plot(x1, y1, color='R', linewidth=1)
    plt.grid(True)
    plt.tight_layout()
    plt.pause(pause_time)
    if islast:
        plt.show()
