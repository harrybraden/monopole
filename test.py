__author__ = 'hwb'

import os
from os import listdir
from os.path import isfile, join
import argparse
import numpy
import copy
from array_tools import smoothed_image, unsmoothed_image
import smoothing_tools
from files import *

from numpy import roots, complex, complex64, mat, dot, trace, pi, sqrt, sum, trace, linalg, matmul, array, matrix, conj, floor, arange
from cmath import exp
import time
from mpmath import ellipk, ellipe, j, taufrom, jtheta, qfrom, ellipf, asin, mfrom
from numpy import roots, complex64, conj, pi, sqrt, sum, trace, linalg, matmul, array
import math
from array_tools import get_point_tuples, load_slice, get_max_value
from array import array

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

directory = "~/Desktop/numerical monopoles/testing_energydensity/"

def write_point_to_file(points, filename):

    """

    :rtype : object
    """
    fo = open(os.path.expanduser(directory + filename), 'wb')
    byteArray = bytearray(points)
    fo.write(byteArray)
    fo.close()

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

    return points

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


def test_higgs_exceptional(k):
    kk = "%1.2f" % k

    points = get_point_tuples(kk, 253, 256)

    exceptional=[]

    for p in points:
        value = higgs_squared(k, p[1], p[0], p[2])
        if (value > 1.0 or value <0.0):
            exceptional.append(p)
            print "Exceptional", p

    write_file = '/Users/hwb/Desktop/numerical monopoles/formaple/higgs_' + str(kk)  +'.txt'

    fo = open(os.path.expanduser(write_file), 'w+')
    for i, point in enumerate(exceptional):
        point_string =  str(point[0]) +' '+ str(point[1]) + ' ' + str(point[2]) +'\n'
        if (i != (len(points) - 1)):
            point_string = point_string
        fo.write(point_string)

    fo.close()

    return

# test_higgs_exceptional(0.90)


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

# print higgs_squared(0.92, 1.175,   0.025, 1.075)


def higgs_on_line(k, x0, x1, y, z, partition_size):

    x_step = (x1 - x0) / partition_size

    for i in range(0, partition_size):
        x=x0+i*x_step
        print higgs_squared(k,x,y,z)

    return


# higgs_on_line(0.8, 0.05, 1.05, 0.1 , 0.5, 50)


def test_on_line(k, x0, x1, y, z, partition_size):

    k1 = sqrt(1-k**2)
    a=k1+complex(0, 1)*k
    b=k1-complex(0, 1)*k

    x_step = (x1 - x0) / partition_size

    for i in range(0, partition_size):
        x=x0+i*x_step
        zeta = calc_zeta(k ,x, y, z)
        eta = calc_eta(k, x, y, z)
        abel = calc_abel(k, zeta, eta)
        mu = calc_mu(k, x, y, z, zeta, abel)

        print "zeta/a", zeta[0]/a, zeta[1]/a, zeta[2]/a, zeta[3]/a, '\n'
        print "eta", eta[0], eta[1], eta[2], eta[3], '\n'
        print "abel",  abel[0], abel[1],abel[2],abel[3],'\n'

        abel_tmp= map(lambda zetai : \
                complex(0, 1) * 1/(complex64(ellipk(k**2))*2*b) \
                * complex64(ellipf( asin( (zetai )/a), mfrom(k=a/b))) \
                - taufrom(k=k)/2,
                zeta)
        print "abel_tmp",  abel_tmp[0], abel_tmp[1],abel_tmp[2],abel_tmp[3],'\n'

        print  abel[0]+conj(abel[2]), abel[1]+conj(abel[3]), '\n'
        print - taufrom(k=k)/2, '\n'
        for l in range(0,4):
            tmp= abs(complex64(calc_eta_by_theta(k, abel[l])) - eta[l])
            print tmp

        value = higgs_squared(k, x, y, z)
        print "higgs", higgs_squared(k,x,y,z)
        if (value > 1.0 or value <0.0):
            print "Exception"

        print '\n'
        # print mu[0]+mu[2], mu[1]+mu[3], '\n'

    return


# test_on_line(0.92, 1.165,  1.195,  0.075, 2.925, 3)

# test_on_line(0.92, 1.175,  1.225,  0.025, 1.075, 1)


# test_on_line(0.92, 0.125, 0.225,  2.625, 1.025,1)

# test_on_line(0.5, 0.5, 1.1, 0.575 , 1.775, 1)

# test_on_line(0.82, 0.9, 1.1, 0.075 , 1.775, 5)

def test_energy_density_exceptional(k):
    kk = "%1.2f" % k

    points = get_point_tuples(kk, 0, 1)

    exceptional=[]

    for p in points:
        value = energy_density(k, p[0], p[1], p[2])
        if (value > 253 ):
            exceptional.append(p)
            print "Exceptional", p, value

    write_file = '/Users/hwb/Desktop/numerical monopoles/formaple/energy_density_' + str(kk)  +'.txt'

    fo = open(os.path.expanduser(write_file), 'w+')
    for i, point in enumerate(exceptional):
        point_string =  str(point[0]) +' '+ str(point[1]) + ' ' + str(point[2]) +'\n'
        if (i != (len(points) - 1)):
            point_string = point_string
        fo.write(point_string)

    fo.close()

    return

# test_energy_density_exceptional(0.95)


# print  len(get_point_tuples(0.80, 253, 256))
# print get_point_tuples(0.80, 253, 256)
# print calc_zeta(0.8, 0.975,0.025, 1.575)
# print energy_density(0.8, 0.975,0.025, 1.575)

def sum_energy_density(k):
    kk = "%1.2f" % k

    directory = '/Users/hwb/Desktop/numerical monopoles/python_results'

    k_directory = directory + '/k=' + kk +'/'

    files = [f for f in listdir(k_directory) if isfile(join(k_directory, f))]

    sum = 0

    for f in files:
        shaped_list = load_slice(k_directory + f)
        for i, xl in enumerate(shaped_list):
            for j, e in enumerate(xl):
                sum=sum+e

    return float(sum)/float(20*20*20*256)


#
# for j in range(1,100):
#     print  sum_energy_density(float(j)/float(100) )


# for j in range(1,20):
#     k1 = float(j)/100.0
#     kk = "%1.2f" % k1
#     print kk
#
#
# print arange(0.95,0.96)



# write_point_to_file(energy_density_on_xz_plane(0.75, 0.025, 3.025, 0.025, 3.025, 0.025, 60), './python_results/floating_test_xz.bin')

#
# points = read_floats('./python_results/floating_test.bin', 25)
# for p in points:
#     print p


# print energy_density(.8, 0.85, 0.0, 0.0)

# print sqrt(higgs_squared(0.001, 0.785, 0.05, 0.05))

# print higgs_squared(0.8 , 6 , 7 , 0.0)

# for i in range(0, 60):
#     x=0.025+i*0.05
#     print x, higgs_squared(0.8,x,1.50,1.50)

# t0 = time.time()
# p = higgs_squared_on_xy_plane(0.8, 0.025, 3.025, 0.025, 3.025, 1.0, 60)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
# t1 = time.time()
# print str(t1-t0)
#
#
# write_point_to_file(p , 'xyhiggs_1.0')


# t0 = time.time()
# p = energy_density_on_xy_plane(0.8, 0.025, 3.025, 0.025, 3.025, 0.1, 60)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
# t1 = time.time()
# print str(t1-t0)
#
#
# write_point_to_file(p , 'xy_0.10')

# print energy_density(0.8, 0.01, 0.0, 0.0)
# #
# print energy_density_at_origin(0.8)
# #
# print energy_density(0.8, 0.3, 0.4, 0.0)

def testing_components_energy_density(k, x1, x2, x3):


    if (is_awc_multiple_root(k, x1, x2, x3) ):
        return 4

    if (is_awc_branch_point(k, x1, x2, x3) ):
        return 4

    zeta = calc_zeta(k ,x1, x2, x3)
    eta = calc_eta(k, x1, x2, x3)
    abel = calc_abel(k, zeta, eta)
    mu = calc_mu(k, x1, x2, x3, zeta, abel)
    x=[x1,x2,x3]

    K = complex64(ellipk(k**2))

    E = complex64(ellipe(k**2))

    cm= (2*E-K)/K

    k1 = sqrt(1-k**2)

    xp = x[0]+complex(0,1)*x[1]
    xm = x[0]-complex(0,1)*x[1]
    S =  sqrt(K**2-4*xp*xm)
    SP = sqrt(K**2-4*xp**2)
    SM = sqrt(K**2-4*xm**2)
    SPM = sqrt(-k1**2*(K**2*k**2-4*xm*xp)+(xm-xp)**2)
    R = 2*K**2*k1**2-S**2-8*x[2]**2
    RM = complex(0,1)*SM**2*(xm*(2*k1**2-1)+xp)-(16*complex(0,1))*xm*x[2]**2
    RP = complex(0,1)*SM**2*(xp*(2*k1**2-1)+xm)+(16*complex(0,1))*xp*x[2]**2
    RMBAR=-complex(0,1)*SP**2*( xp*(2*k1**2-1)+xm ) +16*complex(0,1)*xp*x[2]**2
    RPBAR=-complex(0,1)*SP**2*( xm*(2*k1**2-1)+xp ) -16*complex(0,1)*xm*x[2]**2
    r=sqrt(x[0]**2+x[1]**2+x[2]**2)

    DM = dmus(zeta, x, k)
    DZ = dzetas(zeta, x,k)
    DDM = ddmus(zeta, x, k)
    DDZ = ddzetas(zeta, x,k)

    GNUM = grams(zeta, mu, [x1, x2, x3], k)

    # DGS1 = dgrams1(zeta, mu, DM, DZ, x, k)
    #
    # DGS2 = dgrams2(zeta, mu, DM, DZ, x, k)
    #
    # DGS3 = dgrams3(zeta, mu, DM, DZ, x, k)

    # return  exp(-6 * mu[0])
    print zeta, '\n', abel,  '\n', mu, '\n' #mu[0], abel[0], exp(-6*mu[0])
    return

# print testing_components_energy_density(0.8, 0.975, 0.025, 1.475)
#
# # print testing_components_energy_density(0.8, 0.975, 0.025, 1.540)
# print testing_components_energy_density(0.8, 0.975, 0.025, 1.544)
#
# print testing_components_energy_density(0.8, 0.975, 0.025, 1.575)
#
#
# print energy_density(0.8, 0.975, 0.025, 1.475)
# print energy_density(0.8, 0.975, 0.025, 1.544)
# print int(floor(256*energy_density(0.8, 0.975, 0.025, 1.575)))

# k = '%1.2f' % 0.8
# print k
# directory = '/Users/hwb/Desktop/numerical monopoles/testing_energydensity'
#
# k_directory = directory + '/k=' + str(k) + '/xy_' + '%s' %k + 'y'
#
# print k_directory

# print get_point_tuples(0.8,253,256)

# p= get_point_tuples(0.01,0,0)
#
# print len(p)

# print energy_density(0.99, 1.725, 0.025, 0.025)
#
# for i in range(0, 20):
    # if (energy_density(0.1, p[i][0], p[i][1],p[i][2])==-3):
    # xval= 0.975+0.05*i
    # print  xval, energy_density(0.9, xval, 0.025, 0.025)

# print energy_density_at_origin(0.01)

file = "~/Desktop/numerical monopoles/testing_energydensity/k=0.01/xy_0.01_0.225"

fo   = open(os.path.expanduser(file) )

data = numpy.fromfile(fo)

print len(data)

print data[0]