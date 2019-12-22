from matplotlib import pyplot as plt
import numpy as np


def main(fig, axs, nrow, ncol, x1, y1, x2, y2, pause_time, islast):
    if nrow == 1 and ncol == 1:
        # plt.cla()
        axs.set_xlabel('angle')
        axs.set_ylabel('value')
        axs.set_title('Sin')
        # axs[0, 0].axis('equal')
        axs.set_xlim([0, 10])
        axs.set_ylim([-1.2, 1.2])
        axs.plot(x1, y1, color='B', label='plot 1', linewidth=1)

    elif (nrow == 1 and ncol == 2) or (nrow == 2 and ncol == 1):
        # plt.cla()
        axs[0].set_xlabel('angle')
        axs[0].set_ylabel('value')
        axs[0].set_title('Sin')
        # axs[0, 0].axis('equal')
        axs[0].set_xlim([0, 10])
        axs[0].set_ylim([-1.2, 1.2])
        axs[0].plot(x1, y1, color='B', label='plot 1', linewidth=1)

        # plt.cla()
        axs[1].set_xlabel('angle')
        axs[1].set_ylabel('value')
        axs[1].set_title('Sin')
        # axs[0, 0].axis('equal')
        axs[1].set_xlim([0, 10])
        axs[1].set_ylim([-1.2, 1.2])
        axs[1].plot(x1, y1, color='B', label='plot 1', linewidth=1)
    elif nrow == 2 and ncol == 2:
        # plt.cla()
        axs[0, 0].set_xlabel('angle')
        axs[0, 0].set_ylabel('value')
        axs[0, 0].set_title('Sin')
        # axs[0, 0].axis('equal')
        axs[0, 0].set_xlim([0, 10])
        axs[0, 0].set_ylim([-1.2, 1.2])
        axs[0, 0].plot(x1, y1, color='B', label='plot 1', linewidth=1)

        # plt.cla()
        axs[0, 1].set_xlabel('angle')
        axs[0, 1].set_ylabel('value')
        axs[0, 1].set_title('Cos')
        # axs[0, 1].axis('equal')
        axs[0, 1].set_xlim([0, 10])
        axs[0, 1].set_ylim([-1.2, 1.2])
        axs[0, 1].plot(x2, y2, color='R', label='plot 1', linewidth=1)

        # plt.cla()
        axs[1, 0].set_xlabel('angle')
        axs[1, 0].set_ylabel('value')
        axs[1, 0].set_title('Combine')
        # axs[1, 0].axis('equal')
        axs[1, 0].set_xlim([0, 10])
        axs[1, 0].set_ylim([-1.2, 1.2])
        axs[1, 0].plot(x1, y1, color='B', label='plot 1', linewidth=1)
        axs[1, 0].plot(x2, y2, color='R', label='plot 1', linewidth=1)

        # plt.cla()
        axs[1, 1].set_xlabel('angle')
        axs[1, 1].set_ylabel('difference')
        axs[1, 1].set_title('Differ')
        # axs[1, 1].axis('equal')
        axs[1, 1].set_xlim([0, 10])
        axs[1, 1].set_ylim([-1.2, 1.2])
        axs[1, 1].plot(x1, np.array(y1) - np.array(y2), color='B', label='plot 1', linewidth=1)

    plt.tight_layout()
    plt.pause(pause_time)

    # fig.show()
    if islast:
        plt.show()
    # plt.close(fig)


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
