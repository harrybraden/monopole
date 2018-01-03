__author__ = 'hwb'

import os
from math import *
from numpy import *
from mpmath import *
from os import listdir
from os.path import isfile, join
import time

from higgs_squared import higgs_squared, write_point_to_file, calc_zeta, calc_eta, calc_abel, calc_eta_by_theta, calc_mu

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

def xhiggs_on_line(k, x0, x1, y, z, partition_size):

    x_step = (x1 - x0) / partition_size

    points = []

    for i in range(0, partition_size):
        x=x0+i*x_step
        value = higgs_squared(k, x, y, z)
        bucket_value = int(floor(maxint * value))
        points.append(bucket_value)
        # print higgs_squared(k,x,y,z)

    return points


def yhiggs_on_line(k, x0, x1, y, z, partition_size):

    x_step = (x1 - x0) / partition_size

    points = []

    for i in range(0, partition_size):
        x=x0+i*x_step
        value = higgs_squared(k, y, x, z)
        bucket_value = int(floor(maxint * value))
        points.append(bucket_value)
        # print higgs_squared(k,x,y,z)

    return points


t0 = time.time()
p = higgs_squared_on_xy_plane(0.8, 0.025, 3.025, 0.025, 3.025, 0.05, 60)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
t1 = time.time()
print str(t1-t0)


write_point_to_file(p , 'xyhiggs_0.05')



# px  = xhiggs_on_line(0.8, 0.025, 2.025, 0.10 , 0.025, 50)
#
# py  = yhiggs_on_line(0.8, 0.025, 2.025, 0.10 , 0.025, 50)
# write_point_to_file(p , 'xhiggs_0.025')

# print px
#
# print py

# for i in range(0, 60):
#      print  int(floor(maxint * higgs_squared(0.8, 0.025+0.05*i, 0.025+0.05*(1 -1), 0.10)))

def test_on_line(k, x0, x1, y, z, partition_size):

    k1 = sqrt(1-k**2)
    a=k1+complex(0, 1)*k
    b=k1-complex(0, 1)*k

    # x_step = (x1 - x0) / partition_size
    x_step = 0.05

    for i in range(0, partition_size):
        x=x0+i*x_step
        zeta = calc_zeta(k ,x, y, z)
        eta = calc_eta(k, x, y, z)
        abel = calc_abel(k, zeta, eta)
        mu = calc_mu(k, x, y, z, zeta, abel)

        print i,'\n'

        print "zeta", zeta[0], zeta[1], zeta[2], zeta[3], '\n'
        # print "eta", eta[0], eta[1], eta[2], eta[3], '\n'
        # print "mu", mu[0], mu[1], mu[2], mu[3], '\n'
        # print  mu[0]+conj(mu[2]), mu[1]+conj( mu[3]), '\n'
        # print "abel",  abel[0], abel[1],abel[2],abel[3],'\n'

        # abel_tmp= map(lambda zetai : \
        #                   complex(0, 1) * 1/(complex64(ellipk(k**2))*2*b) \
        #                   * complex64(ellipf( asin( (zetai )/a), mfrom(k=a/b))) \
        #                   - taufrom(k=k)/2,
        #               zeta)
        # print "abel_tmp",  abel_tmp[0], abel_tmp[1],abel_tmp[2],abel_tmp[3],'\n'
        #
        # print  abel[0]+conj(abel[2]), abel[1]+conj(abel[3]), '\n'
        # print - taufrom(k=k)/2, '\n'
        # for l in range(0,4):
        #     tmp= abs(complex64(calc_eta_by_theta(k, abel[l])) - eta[l])
        #     print tmp

        value = higgs_squared(k, x, y, z)
        print "higgs", value, int(floor(maxint *value))
        # if (value > 1.0 or value <0.0):
        #     print "Exception"
        #
        print '\n'
        # # print mu[0]+mu[2], mu[1]+mu[3], '\n'

    return

test_on_line(0.80, 0.025,  1.025,  0.025+0.05*(1-1), 0.1, 2)