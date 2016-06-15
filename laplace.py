from numpy import *
from mpmath import *

from scipy.signal import convolve

z = [0,0,0,0,0]
o = [0,0,-1, 0,0]
s = [0,0,16,0,0]
a = [-1, 16, -90, 16, -1]

stencil = array([[z, z, o, z, z], [z,z,s,z,z], [o,s,a,s,o],[z,z,s,z,z], [z,z,o,z,z]])/12.0


def five_point_laplace(phi, step_size):
    l = convolve(phi, stencil, mode='valid') * float((1.0/step_size)**2)

    return l

#
#
# sinx = []
# for i in range(2, 200, 2):
#     ys = []
#     for j in range(5):
#         zs = []
#         for k in range(5):
#             zs.append(math.sin( float(i)/100 ))
#         ys.append(zs)
#     sinx.append(ys)
#
#
# laplacian = five_point_laplace(sinx, 0.02)
#
#
# print map(lambda si:si[0][0], sinx)
# print laplacian