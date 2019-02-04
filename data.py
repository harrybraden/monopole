"""
Utils for loading data
"""
from numpy import sqrt, int, arange
import numpy as np
import os
import numpy
import copy
import smoothing_tools
import scipy.ndimage

def reconstruct_2d(data):
    n = int(sqrt(data.size))
    points = []
    for i in range(0, n):
        points.append([])
    for i in range(0, data.size):
        y = i / n
        points[y].append(data[i])
    return points
 
def reflect_symmetries(positive_quadrant):
    bottom_right_quadrant = positive_quadrant
    top_right_quadrant = copy.deepcopy(bottom_right_quadrant)
    top_right_quadrant.reverse()

    right_half = top_right_quadrant + bottom_right_quadrant

    left_half = copy.deepcopy(right_half)
    for a in left_half:
        a.reverse()
    full = []
    for i in range(0, len(right_half)):
        full.append(None)
        full[i] = left_half[i] + right_half[i]

    return full

def unsmoothed_data(pth):
    fo = open(os.path.expanduser(pth), 'rb')
    byte_data = numpy.fromfile(fo, dtype=numpy.float64)
    return reconstruct_2d(byte_data)

def smoothed_data(pth):
    return smoothing_tools.smooth_2d(unsmoothed_data(pth))

def reflected_data(pth):
    return reflect_symmetries(unsmoothed_data(pth))#smoothed_data(pth))

def interpolated_data(data, smoothness=3):
    return scipy.ndimage.zoom(data, smoothness, order=3)

def load_data(k, smoothness=3):
    zrange = arange(0.025, 2.975, 0.05)
    points = []
    volume = np.ndarray(shape=(len(zrange)*2,120,120), dtype=float)
    for z, zval in enumerate(zrange):
        pth = './python_smoothed/k=%.02f/xy_%.02f_%.03f' % (k, k, zval)
        dat = reflected_data(pth)
        for y, row in enumerate(dat):
            for x, val in enumerate(row):
                volume[z + len(zrange), y, x] = val
                volume[len(zrange) - z, y, x] = val

    volume = interpolated_data(volume, smoothness=smoothness)
    return volume

