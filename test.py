__author__ = 'hwb'


import os

from numpy import roots, complex, complex64, mat, dot, trace, pi, sqrt, sum, trace, linalg, matmul, array, matrix, conj, floor
from cmath import exp
import time
from mpmath import ellipk, ellipe, j, taufrom, jtheta, qfrom, ellipf, asin, mfrom
from numpy import roots, complex64, conj, pi, sqrt, sum, trace, linalg, matmul, array
import time
import math

from energy_density_old import calc_zeta, calc_eta, calc_abel, calc_mu, energy_density_old, order_roots, quartic_roots
from python_expressions.dmus import dmus
from python_expressions.dzetas import dzetas
from python_expressions.ddmus import ddmus
from python_expressions.ddzetas import ddzetas
from energy_density import energy_density, calc_zeta, calc_eta, calc_abel, calc_mu, order_roots, quartic_roots

# The files above are common in a calculation. They are calculated once and used numerous times.
# The file below is essentially the very long line and whose speed is the question.
from python_expressions.ddphis111 import ddphis111
from python_expressions.NDD111 import nddphis111
from python_expressions.ddgrams211 import ddgrams211


def energy_density_on_xy_plane(k, x0, x1, y0, y1, z, partition_size):

    x_step = (x1 - x0) / partition_size
    y_step = (y1 - y0) / partition_size

    points = []
    last = 0

    for j in range(0, partition_size):
        for i in range(0, partition_size):
            x = x0 + i * x_step
            y = y0 + j * y_step

            value = energy_density(k, x, y, z)
            bucket_value = int(floor(256 * value))
            if(bucket_value > 255 or bucket_value < 0):
                print i, j, bucket_value
                bucket_value = 255
            points.append(bucket_value)

            last = bucket_value

    return points

def energy_density_on_yz_plane(k, y0, y1, z0, z1, x, partition_size):

    y_step = (y1 - y0) / partition_size
    z_step = (z1 - z0) / partition_size

    points = []
    last = 0
    for j in range(0, partition_size):
        for i in range(0, partition_size):
            y = y0 + i * y_step
            z = z0 + j * z_step

            value = energy_density(k, x, y, z)
            bucket_value = int(floor(256 * value))
            if(bucket_value > 255 or bucket_value < 0):
                print i, j, bucket_value
                bucket_value = 255
            points.append(bucket_value)

            last = bucket_value

    return points

def energy_density_on_xz_plane(k, x0, x1, z0, z1, y, partition_size):

    x_step = (x1 - x0) / partition_size
    z_step = (z1 - z0) / partition_size


    points = []
    last = 0
    for j in range(0, partition_size):
        for i in range(0, partition_size):
            x = x0 + i * x_step
            z = z0 + j * z_step

            value = energy_density(k, x, y, z)
            bucket_value = int(floor(256 * value))
            if(bucket_value > 255 or bucket_value < 0):
                print i, j, bucket_value
                bucket_value = 255
            points.append(bucket_value)

            last = bucket_value

    return points


def write_point_to_file(points, filename):

    """

    :rtype : object
    """
    fo = open(os.path.expanduser("~/Desktop/numerical monopoles/python_results/" + filename), 'wb')
    byteArray = bytearray(points)
    fo.write(byteArray)
    fo.close()

# t1 = time.time()
# p = energy_density_on_xy_plane(0.999, 0.1, 3.1, 0.1, 3.1, 0, 30)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
# t2 = time.time()
#
# print str(t2-t1)
#

# for i in range(0, 5):
#     z = float(i) / 2
#     p = energy_density_on_xz_plane(0.8, 0.1, 3.1, 0.1, 3.1, z, 40)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
#     write_point_to_file(p , 'example_xz_8_3_' + str(i))
#     print i

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
#




# print energy_density(0.8, 2.6, 1.9, 0)  # I find the energy density here appears a minimum
#
# print int(floor(256 * energy_density(0.8, 2.6, 1.9, 0)))


def energy_density_at_origin(k):
    K = complex64(ellipk(k**2))
    E = complex64(ellipe(k**2))
    k1 = sqrt(1-k**2)

    A = 32*(k**2 *(-K**2 * k**2 +E**2-4*E*K+3* K**2 + k**2)-2*(E-K)**2)**2/(k**8 * K**4 * k1**2)

    return A.real


# for i in range(0, 10 , 1):
#     k = 0.05+i* 0.1
#     print k
#     print energy_density_at_origin(k)
#     print energy_density(k, 0.01, 0, 0)
#     print energy_density(k, 0, 0.01, 0)
#     print energy_density(k, 0, 0, 0.01)




# p = energy_density_on_xz_plane(0.8, 0.1, 3.1, 0.1, 3.1, 0.0, 30)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
#
# write_point_to_file(p , 'example_xz' + str(0.0))

# print energy_density(0.8, 1, 0.1, 0.1)
# print energy_density(0.8, 1, 3.1, 0.1)
# print energy_density(0.8, 1, 0.1, 3.1)

# for j in range(0, 40):
#     for i in range(0, 40):
#         x = 0.1 + i * 0.075
#         z = 0.1 + j * 0.075
#
#         print i,j, energy_density(0.8, x, 0, z)


# energy_density_on_xz_plane(0.8, 0.1, 3.1, 0.1, 3.1, 0, 40)
#  k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final



def test_timing(k, x1, x2, x3):

    t2 =  time.time()

    A = energy_density_old(k, x1, x2, x3)

    t3= time.time()

    B = energy_density(k, x1, x2, x3)

    t4= time.time()

    print str(t3-t2)
    print str(t4-t3)


    return A, B
#
# t15 = time.time()
# C  = test_timing(.2, 1.5, 1.2, 0.6)
# t16 = time.time()
#
# print C
# print str(t16-t15)



# for j in range(0, 2):
#     for i in range(0, 2):
#         x = 0.1 + i * 0.075
#         y = 0.1 + j * 0.075
#
#         print test_timing(0.8, x, y, 1.0)

# energy_density_on_xy_plane(0.8, 0.1, 2.1, 0.1, 2.1, 0.1, 40)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final

def is_awc_multiple_root(k, x1, x2, x3):
    K = complex64(ellipk(k**2))
    k1 = sqrt(1-k**2)
    tol = 0.001

    if (4 * x1**2 * k1**2 - k**2 *( K**2 *k1**2 + 4* x2**2) < tol) and x3==0:
        return True
    elif (4 * x1**2 * k1**2 - K**2 *k1**2 + 4* x3**2 < tol) and x2==0:
        return True

    return False



# print energy_density(0.8, 1.0, 0, 2.35)

# p = energy_density_on_xy_plane(0.8, 0.1, 2.1, 0.1, 2.1, 0.0, 20)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final

print energy_density(0.8, 0.4, 0.3, 0.0)

print is_awc_multiple_root(0.8, 0.4, 0.3, 0.0)