from matplotlib import pyplot as plt
import numpy as np


def main(x1, y1, x2, y2, pause_time, islast):
    plt.figure(1)
    plt.clf()
    plt.xlabel('age')
    plt.ylabel('salary')
    plt.title('Sin')
    plt.plot(x1, y1, label='plot 1', linewidth=1)
    plt.grid(False)
    plt.tight_layout()
    plt.pause(pause_time)
    if islast:
        plt.show()

    plt.figure(2)
    plt.clf()
    plt.xlabel('age')
    plt.ylabel('salary')
    plt.title('Cos')
    plt.plot(x2, y2, label='plot 2', linewidth=1)
    plt.grid(False)
    plt.tight_layout()
    plt.pause(pause_time)
    if islast:
        plt.show()
iterate = 10
a = list(np.arange(0,iterate, 1))
print(a)
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