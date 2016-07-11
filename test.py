__author__ = 'hwb'

import time
from math import *
from mpmath import *
import os

from energy_density import calc_zeta, calc_eta, calc_abel, calc_mu, energy_density
from python_expressions.dmus import dmus
from python_expressions.dzetas import dzetas
from python_expressions.ddmus import ddmus
from python_expressions.ddzetas import ddzetas

# The files above are common in a calculation. They are calculated once and used numerous times.
# The file below is essentially the very long line and whose speed is the question.
from python_expressions.ddphis111 import ddphis111
from python_expressions.NDD111 import nddphis111
from python_expressions.ddgrams211 import ddgrams211

def test_timing(k, x1, x2, x3):

    t0 = time.time()
    zeta = calc_zeta(k ,x1, x2, x3)
    eta = calc_eta(k, x1, x2, x3)
    abel = calc_abel(k, zeta, eta)
    mu = calc_mu(k, x1, x2, x3, zeta, abel)
    x=[x1,x2,x3]


    t1 = time.time()



    DM = dmus(zeta, x, k)
    DZ = dzetas(zeta, x,k)
    DDM = ddmus(zeta, x, k)
    DDZ = ddzetas(zeta, x,k)

    t2 =  time.time()

    A = ddphis111(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

    t3= time.time()

    B = nddphis111(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

    t4= time.time()

    C = ddgrams211(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

    t5= time.time()


    print str(t1-t0)
    print str(t2-t1)
    print str(t3-t2)
    print str(t4-t3)
    print str(t5-t4)



    return A, B # ddphis111(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

# t15 = time.time()
# C  = test_timing(.8, 1.5, 0.5, 0.2)
# t16 = time.time()
#
# print C
# print str(t16-t15)


# print energy_density(0.8, 1/float(2*3) , 0.0 , 0.0)


def energy_density_on_a_plane(k, x0, x1, y0, y1, z, partition_size):

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
                bucket_value = last
            points.append(bucket_value)

            last = bucket_value

    return points



def write_point_to_file(points, filename):

    fo = open(os.path.expanduser("~/Desktop/numerical monopoles/" + filename), 'wb')
    byteArray = bytearray(points)
    fo.write(byteArray)
    fo.close()

t1 = time.time()
p = energy_density_on_a_plane(0.8, 0.1, 5.1, 0.1, 5.1, 0, 100)   # k, x-initial x-final, y-initial, y-final, z, partition size=no points between initial final
t2 = time.time()

print str(t2-t1)

write_point_to_file(p , 'small5')

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
