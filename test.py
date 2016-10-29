__author__ = 'hwb'

import os
import numpy
import copy
from array_tools import smoothed_image, unsmoothed_image
import smoothing_tools

from numpy import roots, complex, complex64, mat, dot, trace, pi, sqrt, sum, trace, linalg, matmul, array, matrix, conj, floor
from cmath import exp
import time
from mpmath import ellipk, ellipe, j, taufrom, jtheta, qfrom, ellipf, asin, mfrom
from numpy import roots, complex64, conj, pi, sqrt, sum, trace, linalg, matmul, array
import math

from higgs_squared import higgs_squared

# from energy_density_old import calc_zeta, calc_eta, calc_abel, calc_mu, energy_density_old, order_roots, quartic_roots
from python_expressions.dmus import dmus
from python_expressions.dzetas import dzetas
from python_expressions.ddmus import ddmus
from python_expressions.ddzetas import ddzetas
from energy_density import energy_density, calc_zeta, calc_eta, calc_abel, calc_mu, order_roots, quartic_roots, \
     energy_density_on_xy_plane, energy_density_on_yz_plane, energy_density_on_xz_plane, energy_density_at_origin, is_awc_multiple_root, \
     is_awc_branch_point, write_point_to_file
from energy_density_old import energy_density_old

# Some global variables.

directory = "~/Desktop/numerical monopoles/python_results/"


# The files above are common in a calculation. They are calculated once and used numerous times.
# The file below is essentially the very long line and whose speed is the question.
from python_expressions.ddphis111 import ddphis111
from python_expressions.NDD111 import nddphis111
from python_expressions.ddgrams211 import ddgrams211


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

# print energy_density(0.8, 1.0, 0, 2.35)

# p = energy_density_on_xy_plane(0.8, 0.1, 3.1, 0.1, 3.1, 0.1, 40)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
#
# write_point_to_file(p, 'example_xy_8_3_1')
#
# for i in range(0, 17):
#     print i
#     z = float(i)*0.1
#     p = energy_density_on_xy_plane(0.2, 0.1, 3.1, 0.1, 3.1, z, 40)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
#     write_point_to_file(p , 'example_xy_2_3_' +str(i))

A = time.time()

p = energy_density_on_xy_plane(0.75, 0.1, 3.1, 0.1, 3.1, 0.2, 10)   # k, x-initial z-final, z-initial, y-final, y, partition size=no points between initial final
filename= (directory + 'example_xy_75_3_2')
write_point_to_file(p ,filename )

B =time.time()
print str(B-A)


# for i in range(0, 17):
#     print i
#     y = float(i)*0.1
#     p = energy_density_on_xz_plane(0.2, 0.1, 3.1, 0.1, 3.1, y, 40)   # k, x-initial z-final, z-initial, y-final, y, partition size=no points between initial final
#     write_point_to_file(p , 'example_xz_2_3_' +str(i))

# p = energy_density_on_xz_plane(0.8, 0.05, 1.05, 0.05, 1.05, 0.0, 10)   # k, x-initial z-final, z-initial, y-final, y, partition size=no points between initial final
# write_point_to_file(p , 'example_xz_8_1_0')

# for i in range(0, 17):
#     print i
#     x = float(i)*0.1
#     p = energy_density_on_yz_plane(0.2, 0.1, 3.1, 0.1, 3.1, x, 40)   # k, y-initial y-final, z-initial, z-final, x, partition size=no points between initial final
#     write_point_to_file(p , 'example_yz_2_3_' +str(i))



# print 'example_xy_8_3_' +str(1) +'v'

# write_point_to_file(energy_density_on_xy_plane(0.8, 0.1, 3.1, 0.1, 3.1, 0.0, 40),  'example_xy_8_3_0v')

# energy_density_on_xy_plane(0.8, 0.1, 2.1, 0.1, 2.1, 0.05, 20)

# print order_roots(quartic_roots(0.8, 0.4, 0.3, 0.0))
# print energy_density(0.8, 0.4, 0.3, 0.0)
# print calc_eta(0.8, 0.4, 0.3, 0.0)



# print order_roots(quartic_roots(0.8, 0.8, 0.6, 0.0))
# print order_roots(quartic_roots(0.8, 1.2, 0.9, 0.0))
# print order_roots(quartic_roots(0.8, 1.6, 1.2, 0.0))
# print order_roots(quartic_roots(0.8, 2.0, 1.5, 0.0))

# p= energy_density_on_xz_plane(0.2, 0.1, 3.1, 0.1, 3.1, 0.0, 40)
#
# write_point_to_file(p, 'example_xz_2_3_0')

# for i in range(5, 8):
#     print i
#     y = float(i)*0.1
#     p = energy_density_on_xz_plane(0.2, 0.1, 3.1, 0.1, 3.1, y, 40)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
#     write_point_to_file(p , 'example_xz_2_3_' + str(i))
#


# k  = 0.8
# x1 = 1.0
# x2 = 0.00
# x3 = 0.0
#
# k1 = sqrt(1-k**2)
# a=k1+complex(0, 1)*k
# b=k1-complex(0, 1)*k


# zeta = calc_zeta(k ,x1, x2, x3)
# eta = calc_eta(k, x1, x2, x3)
# abel = calc_abel(k, zeta, eta)
# mu = calc_mu(k, x1, x2, x3, zeta, abel)
#
# lambda_test = map(lambda zetai : \
#                     complex64(ellipf( asin( (zetai )/a), mfrom(k=a/b))), \
#                 zeta)
# print zeta
# print lambda_test

# for i in range(0, 7, 1):
#     x1 = 0.97+float(i)*0.01
#     zeta = calc_zeta(k ,x1, x2, x3)
#     # print x1, map(lambda zetai :  complex(0, 1) * 1/(complex64(ellipk(k**2))*2*b) *\
#     #                              complex64(ellipf( asin( (zetai )/a), mfrom(k=a/b))),\
#     #             zeta)
#     print x1, zeta[0], map(lambda zetai :  complex(0, 1) * 1/(complex64(ellipk(k**2))*2*b) * \
#                            complex64(ellipf( asin( (zetai )/a), mfrom(k=a/b))),zeta)[0]


# print zeta
# print eta
# print abel
# print mu
#
#
# print higgs_squared(k,x1,x2,x3)
#
# print energy_density_old(k,x1,x2,x3)
#
# print energy_density(k,x1,x2,x3)


# print is_awc_multiple_root(0.8, 1.0, 0.0, 0.1)
#
# print is_awc_branch_point(0.8, 1.0, 0.0, 0.1)

# print energy_density(0.8, 0.4, 0.3, 0.0)
#
# print higgs_squared(0.8, 0.4, 0.3, 0.0)

# data = "example_xy_8_3_1"
#
# file = "~/Desktop/numerical monopoles/python_results/%s"  % data
#
# print file

# def naming(file):
#     data ="file"
#     return "~/Desktop/numerical monopoles/python_results/%s"  % file
#
# print naming("example_xy_8_3_1")


