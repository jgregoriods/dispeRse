import matplotlib.pyplot as plt
import numpy as np

from ctypes import *


lib_abm = CDLL('abm.so')


def run_model(height, width, grid, num_iter):
    c_height = (POINTER(c_int))(c_int(height))
    c_width = (POINTER(c_int))(c_int(width))
    c_grid = (c_int * len(grid))(*grid)
    c_num_iter = (POINTER(c_int))(c_int(num_iter))
    lib_abm.run_model(c_height, c_width, c_grid, c_num_iter)
    return np.ctypeslib.as_array(c_grid)


if __name__ == '__main__':
    grd = [0] * 25 * 50
    res = run_model(50, 25, grd, 500)
    plt.imshow(res.reshape((50, 25)))
    plt.show()
