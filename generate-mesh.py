from numpy import sqrt, int, arange
import numpy as np
import os
import numpy
import copy
import smoothing_tools
from scipy.spatial import Delaunay, ConvexHull

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
    byte_data = numpy.fromfile(fo, dtype=numpy.uint8)
    return reconstruct_2d(byte_data)

def smoothed_data(pth):
    return smoothing_tools.smooth_2d(unsmoothed_data(pth))

def reflected_data(pth):
    return reflect_symmetries(smoothed_data(pth))

def load_data(k, intensity):
    zrange = arange(0.025, 2.975, 0.05)
    points = []
    for z, zval in enumerate(zrange):
        pth = './python_results/k=%.02f/xy_%s_0.025-3.025_0.025-3.025_%s_60' % (k, k, zval)
        dat = reflected_data(pth)
        for y, row in enumerate(dat):
            ''' Scanline through the row, storing entry and exit points'''
            inside = False
            for x, val in enumerate(row):
                if (not inside and val >= intensity) or (inside and val < intensity):
                    points.append([x,y, z])
                    points.append([x,y, -z]) # Symmetrical through Z
                    inside = not inside

    npoints = np.array(points)
    tri = Delaunay(points, qhull_options='Qv Qz')
    #print tri.simplices[0]
    print "# OBJ file"

    for v in tri.points:
        print ("v %d %d %d" % (v[0], v[1], v[2]))

    for f in tri.simplices:
        print ("f %d %d %d" % (f[1] + 1, f[2] + 1, f[3] + 1))

load_data(0.4, 150)
