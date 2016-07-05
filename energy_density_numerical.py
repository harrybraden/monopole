__author__ = 'hwb'
#  This will calculate using a numerical laplacian the energy density based on the higgs_squared
#
from numpy import array
from mpmath import *
import time
import matrices
import os
# from laplace import five_point_laplace
import math
from higgs_squared import higgs_squared

from scipy.signal import convolve




def five_point_laplace(phi, step_size):

    z = [0,0,0,0,0]
    o = [0,0,-1, 0,0]
    s = [0,0,16,0,0]
    a = [-1, 16, -90, 16, -1]

    stencil = array([[z, z, o, z, z], [z,z,s,z,z], [o,s,a,s,o],[z,z,s,z,z], [z,z,o,z,z]])/12.0

    l = convolve(phi, stencil, mode='valid') * float((1.0/step_size)**2)

    return l

def energy_density_numerical(k, x1, x2, x3):
    step_size = 0.02

    points = []
    for a in range(-2, 3, 1):
        points_y = []
        for b in range(-2, 3, 1):
            points_z = []
            for c in range(-2, 3, 1):
                points_z.append(higgs_squared(k, float(x1 + a * step_size), float(x2 + b * step_size), float(x3 + c * step_size)))
            points_y.append(points_z)
        points.append(points_y)

    return five_point_laplace(points, step_size)[0][0]


print energy_density_numerical(0.8, 1.5 , 0.00, 0.0)[0]

    # for i in range(0, 20, 1):
    #     print energy_density_numerical(0.8, (float(2*i) +90)/100 , 0, 0)

def energy_density_on_line(k, x0, y0, z0, axis, end):
    if (axis not in ['x','y','z']):
        raise ValueError("Invalid axis given. Must be one of 'x', 'y', or 'z'")

    step_size = 0.02
    if (axis == 'x'):
        start = x0
    elif (axis == 'y'):
        start = y0
    elif (axis == 'z'):
        start = z0
    intervals = int(math.ceil((end - start)/step_size))


    points = []
    for a in range(-2, intervals + 2, 1):
        points_y = []
        for b in range(-2, 3, 1):
            points_z = []
            for c in range(-2, 3, 1):
                if(b == 0 or c == 0):
                    if (axis == 'x'):
                        points_z.append(higgs_squared(k, float(x0 + a * step_size), float(y0 + b * step_size), float(z0 + c * step_size)))
                    elif (axis == 'y'):
                        points_z.append(higgs_squared(k, float(x0 + b * step_size), float(y0 + a * step_size), float(z0 + c * step_size)))
                    elif (axis == 'z'):
                        points_z.append(higgs_squared(k, float(x0 + c * step_size), float(y0 + b * step_size), float(z0 + a * step_size)))
                else:
                    points_z.append(0) # Value is not used in laplace calculation
            points_y.append(points_z)
        points.append(points_y)

    return five_point_laplace(points, step_size)[0][0]

def energy_density_volume(k, x0, x1, y0, y1, z0, z1):
    step_size = 0.02

    xintervals = int(math.ceil((x1 - x0)/step_size))
    yintervals = int(math.ceil((y1 - y0)/step_size))
    zintervals = int(math.ceil((z1 - z0)/step_size))

    points = []
    for a in range(-2, xintervals + 2, 1):
        points_y = []
        for b in range(-2, yintervals + 2, 1):
            points_z = []
            for c in range(-2, zintervals + 2, 1):
                points_z.append(higgs_squared(k, float(x0 + a * step_size), float(y0 + b * step_size), float(z0 + c * step_size)))
            points_y.append(points_z)
        points.append(points_y)

    return five_point_laplace(points, step_size)



# fo = open(os.path.expanduser("~/Desktop/numerical monopoles/hwb_zlaplaceR"), 'w' )
# ed = energy_density_on_line(0.8, 0, 0, .01, 'z', 4.01)
# for i in range(0, len(ed), 1):
#     fo.write("%4.3f %15.9f\n"% ( 0.01 + (4.01-0.01)*(float(i)/(len(ed)-1)), ed[i]))
# fo.close()

# fo = open(os.path.expanduser("~/Desktop/hwb_xlaplaceR"), 'w' )
# ed = energy_density_on_line(0.8, .005, 0, 0, 'x', 4.005)
# for i in range(0, len(ed), 1):
#     fo.write("%4.3f %15.9f\n"% ( 0.005 + (4.005-0.005)*(float(i)/(len(ed)-1)), ed[i]))
# fo.close()

# print energy_density_on_line(0.8, 4, 0, 0, 'x', 5)

# print energy_density_volume(0.8, 0.5, 0.6, 0.5, 0.6, 0.5, 0.6)
