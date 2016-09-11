from numpy import sqrt, int
import os
import numpy
import copy
import smoothing_tools

def flatten_2d(data):
    return None


def flatten_3d(data):
    return None


def reconstruct_2d(data):
    n = int(sqrt(data.size))
    points = []
    for i in range(0, n):
        points.append([])
    for i in range(0, data.size):
        y = i / n
        points[y].append(data[i])
    return points

def reconstruct_3d(data):
    return None

def smoothed_image(data):
    file = "~/Desktop/numerical monopoles/python_results/%s"  % data

    fo = open(os.path.expanduser(file), 'rb')

    bytes = numpy.fromfile(fo, dtype=numpy.uint8)

    shaped_list = reconstruct_2d(bytes)

    smoothed = smoothing_tools.smooth_2d(shaped_list)

    bottom_right_quadrant = smoothed
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

def unsmoothed_image(data):
    file = "~/Desktop/numerical monopoles/python_results/%s"  % data

    fo = open(os.path.expanduser(file), 'rb')

    bytes = numpy.fromfile(fo, dtype=numpy.uint8)

    shaped_list = reconstruct_2d(bytes)

    bottom_right_quadrant = shaped_list
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

