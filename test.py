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
from array_tools import get_point_tuples

from higgs_squared import higgs_squared


# from energy_density_old import calc_zeta, calc_eta, calc_abel, calc_mu, energy_density_old, order_roots, quartic_roots
from python_expressions.dmus import dmus
from python_expressions.dzetas import dzetas
from python_expressions.ddmus import ddmus
from python_expressions.ddzetas import ddzetas
from energy_density import energy_density, calc_zeta, calc_eta, calc_abel, calc_mu, calc_eta_by_theta, order_roots, quartic_roots, \
     energy_density_on_xy_plane, energy_density_on_yz_plane, energy_density_on_xz_plane, energy_density_at_origin, is_awc_multiple_root, \
     is_awc_branch_point, write_point_to_file
from energy_density_old import energy_density_old

from python_expressions.grams import grams
from modified_expressions.dgrams1 import dgrams1
from modified_expressions.dgrams2 import dgrams2
from modified_expressions.dgrams3 import dgrams3


# Some global variables.

directory = "~/Desktop/numerical monopoles/exceptional/"


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

# p = energy_density_on_xy_plane(0.75, 0.025, 3.025, 0.025, 3.025, 1.0 , 60)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
#
# write_point_to_file(p, 'example_xy_8_3_1')
#
# for i in range(0, 17):
#     print i
#     z = float(i)*0.1
#     p = energy_density_on_xy_plane(0.2, 0.1, 3.1, 0.1, 3.1, z, 40)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
#     write_point_to_file(p , 'example_xy_2_3_' +str(i))

# A = time.time()
#
# p = energy_density_on_xy_plane(0.80, 0.95, 1.05, 0.025, 3.025, 3.0 , 20)   # k, x-initial z-final, z-initial, y-final, y, partition size=no points between initial final
# filename= (directory + 'example_xy_80')
# write_point_to_file(p ,filename )
#
# B =time.time()
# print str(B-A)


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


k  = 0.8
K = complex64(ellipk(k**2))
x1 = K/2
x2 = 0.02
x3 = 2.0

k1 = sqrt(1-k**2)
a=k1+complex(0, 1)*k
b=k1-complex(0, 1)*k


zeta = calc_zeta(k ,x1, x2, x3)
eta = calc_eta(k, x1, x2, x3)
abel = calc_abel(k, zeta, eta)
mu = calc_mu(k, x1, x2, x3, zeta, abel)
x=[x1,x2,x3]

DM = dmus(zeta, x, k)
DZ = dzetas(zeta, x,k)
DDM = ddmus(zeta, x, k)
DDZ = ddzetas(zeta, x,k)


#
# lambda_test = map(lambda zetai : \
#                     complex64(ellipf( asin( (zetai )/a), mfrom(k=a/b))), \
#                 zeta)
# print zeta
# print lambda_test

# for i in range(1, 7, 1):
#     x2 = float(i)*0.001
#
#     print x2
#
#     print higgs_squared(k,x1,x2,x3)
#
#     print energy_density(k,x1, x2,x3)

    # print x1, map(lambda zetai :  complex(0, 1) * 1/(complex64(ellipk(k**2))*2*b) *\
    #                              complex64(ellipf( asin( (zetai )/a), mfrom(k=a/b))),\
    #             zeta)
    # print x1 , zeta[0], map(lambda zetai :  complex(0, 1) * 1/(complex64(ellipk(k**2))*2*b) * \
              #             complex64(ellipf( asin( (zetai )/a), mfrom(k=a/b))),zeta)[0]

# print x1
# print zeta
# print eta
# print abel
# print mu
# print higgs_squared(k,x1,x2,x3)
# print DM
# print DDM
# print DZ
# print DDZ
# print grams(zeta, mu, [x1, x2, x3], k)
#
# # print phis(zeta, mu, [x1, x2, x3], k)
# print energy_density(k,x1,x2,x3)
# print '\n'
#
# print DZ
# print '\n'
#
# print calc_zeta(k ,x1, 0.0, x3)
# print energy_density(k,x1,0.0,x3)


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


maxint = 256                     # This determines the digits being returned
maxintr = maxint -1


def higgs_squared_on_xy_plane(k, x0, x1, y0, y1, z, partition_size):  # If this falls outside of [0,1) an value of maxint-1 will be returned; maxint is globally set.

    x_step = (x1 - x0) / partition_size
    y_step = (y1 - y0) / partition_size

    points = []

    for j in xrange(0, partition_size):
        # if j % 10 == 0 and j > 0:
        #     eprint("- rendered %s lines..." % j)

        for i in xrange(0, partition_size):
            x = x0 + i * x_step
            y = y0 + j * y_step

            value = higgs_squared(k, x, y, z)
            bucket_value = int(floor(maxint * value))
            if(bucket_value > maxintr or bucket_value < 0):
                print i, j, bucket_value
                bucket_value = maxintr
            points.append(bucket_value)

def higgs_squared_on_yz_plane(k, y0, y1, z0, z1, x, partition_size):   # If this falls outside of [0,1) an value of maxint-1 will be returned; maxint is globally set.

    y_step = (y1 - y0) / partition_size
    z_step = (z1 - z0) / partition_size

    points = []
    for j in range(0, partition_size):
        for i in range(0, partition_size):
            y = y0 + i * y_step
            z = z0 + j * z_step

            value = higgs_squared(k, x, y, z)
            bucket_value = int(floor(maxint * value))
            if(bucket_value > maxintr or bucket_value < 0):
                print i, j, bucket_value
                bucket_value = maxintr
            points.append(bucket_value)


    return points


    return points


def test_higgs_exceptional(k):
    kk = "%1.2f" % k

    points = get_point_tuples(kk, 253, 256)

    exceptional=[]

    for p in points:
        value = higgs_squared(k, p[1], p[0], p[2])
        if (value > 1.0 or value <0.0):
            exceptional.append(p)

    write_file = '/Users/hwb/Desktop/numerical monopoles/formaple/higgs_' + str(kk)  +'.txt'

    fo = open(os.path.expanduser(write_file), 'w+')
    for i, point in enumerate(points):
        point_string =  str(point[0]) +' '+ str(point[1]) + ' ' + str(point[2]) +'\n'
        if (i != (len(points) - 1)):
            point_string = point_string
        fo.write(point_string)

    fo.close()

    return

# test_higgs_exceptional(0.82)


# def test_print(k):
#     kk = "%1.2f" % k
#     print  '/Users/hwb/Desktop/numerical monopoles/formaple/higgs_' + str(kk)  +'.txt'
#
#
# test_print(0.92)

def zeta_on_line(k, x0, x1, y, z, partition_size):

    x_step = (x1 - x0) / partition_size

    for i in range(0, partition_size):
        x=x0+i*x_step
        print calc_zeta(k,x,y,z)

    return

# zeta_on_line(0.82, 0.9, 1.1, 0.075 , 1.775, 50)

# print higgs_squared(0.82, 0.975, 0.075 , 1.775)

def higgs_on_line(k, x0, x1, y, z, partition_size):

    x_step = (x1 - x0) / partition_size

    for i in range(0, partition_size):
        x=x0+i*x_step
        print higgs_squared(k,x,y,z)

    return


# higgs_on_line(0.82, 0.9, 1.1, 0.075 , 1.775, 50)


def test_on_line(k, x0, x1, y, z, partition_size):

    x_step = (x1 - x0) / partition_size

    for i in range(0, partition_size):
        x=x0+i*x_step
        zeta = calc_zeta(k ,x, y, z)
        eta = calc_eta(k, x, y, z)
        abel = calc_abel(k, zeta, eta)
        mu = calc_mu(k, x, y, z, zeta, abel)
        smu =mu[0]+mu[1]+mu[2]+mu[3]
        print mu[0]+mu[2], mu[1]+mu[3], '\n'
        print abel[0]+conj(abel[2]), abel[1]+conj(abel[3]), '\n'
        print - taufrom(k=k)/2
        for k in range(0,3):
            tmp= abs(complex64(calc_eta_by_theta(k, abel[k])) - eta[k])
            print tmp

    return


# test_on_line(0.82, 0.9, 1.1, 0.075 , 1.775, 50)


test_on_line(0.5, 0.5, 1.1, 0.575 , 0.775, 1)