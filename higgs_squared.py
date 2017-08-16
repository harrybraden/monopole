__author__ = 'hwb'
#  This file will calculate the trace of the higgs field squared for a point (x1, x2, x3) of space and a parameter k (between 0 and 1)
#
#  Given a point x we associate to this four points P=(zeta, eta) on an elliptic curve. The four zeta values for these are determined by the roots of the quartic describing the curve
#  quartic_roots(k, x1, x2, x3) but they are ordered by properties coming from the real structure order_roots(roots); then
#  calc_zeta(k, x1, x2, x3)= order_roots(quartic_roots(k, x1, x2, x3))
#  The second coordinate eta is given by calc_eta
#
#  We need to calculate the Abel image of P. This will be done in calc_abel. To correctly identify this point with the correct sheet we will calculate the eta value of the Abel image by
#  calc_eta_by_theta and if this agrees accept it, and if not shift by the half period to the correct sheet. This is done by abel_select and calc_abel returns the correct Abel images
#  for each of the points.
#
#  There is one transcendental function mu for each of the four points P given by calc_mu(k, x1, x2, x3, zeta, abel).
#
#  All functions are then functions of (k, x1, x2, x3) and zeta_i, eta_i, mu_i (i=1..4).
#

from numpy import roots, complex, complex64, complex128, mat, dot, trace, pi, sqrt, sum, trace, linalg, matmul, array, matrix, conj, sort_complex
from math import *
from cmath import exp
import time
from mpmath import ellipk, ellipe, j, taufrom, jtheta, qfrom, ellipf, asin, mfrom
from numpy import roots, complex64, conj, pi, sqrt, sum, trace, linalg, matmul, array
import time
import math
import os





from python_expressions.grams import grams
from python_expressions.phis import phis

def quartic_roots(k, x1, x2, x3):
    K = complex128(ellipk(k**2))

    e0 = complex128((x1*j - x2)**2 + .25 * K**2)
    e1 = complex128(4*(x1*j-x2)*x3)
    e2 = complex128(4*(x3**2) - 2 * (x1**2) - 2 * (x2**2) + (K**2) * (k**2 - 0.5))
    e3 = complex128(4*x3*(x2 + j*x1))
    e4 = complex128(x2**2 - x1**2 + 2*j*x1*x2 + 0.25*K**2)

    return sort_complex(roots([e4, e3, e2, e1, e0]))     # I put the sort_complex to have a canonical form, so that when we order them they will vary continuously
    # return roots([e4, e3, e2, e1, e0])


def order_roots(roots):
    if len(roots) != 4:
        raise ValueError

    if abs( roots[0]+1/conj(roots[1]) ) <10**(-4):
        return [roots[0], roots[2], roots[1], roots[3]]
    elif abs(roots[0]+1/conj(roots[2])) <10**(-4):
        return [roots[0],roots[1],roots[2],roots[3]]
    elif abs(roots[0]+1/conj(roots[3])) <10**(-4):
        return[roots[0],roots[1],roots[3],roots[2]]
    raise ValueError

def calc_zeta(k, x1, x2, x3):
    return order_roots(quartic_roots(k, x1, x2, x3))

def calc_eta(k, x1, x2, x3):
    zeta = calc_zeta(k, x1, x2, x3)
    return map(lambda zetai : -(x2 + j*x1) * zetai**2 - 2*x3 * zetai + x2 - j*x1, zeta)

def calc_abel(k, zeta, eta):
    k1 = sqrt(1-k**2)

    a=k1+complex(0, 1)*k

    b=k1-complex(0, 1)*k

    abel_tmp = map(lambda zetai : \
                       complex(0, 1) * 1/(complex64(ellipk(k**2))*2*b) \
                       * complex64(ellipf( asin( (zetai )/a), mfrom(k=a/b))) \
                       - taufrom(k=k)/2,
                   zeta)

    abel = []
    for i in range(0, 4, 1):
        abel.append(abel_select(k, abel_tmp[i], eta[i]))

    return abel

def abel_select(k, abeli, etai):
    # tol = 0.001
    # tol = 1.0                 # This choice appears quite important to get smoothness. It was too small initially.
    tol = 0.10                 # This choice appears quite important to get smoothness. It was too small initially.


    if (abs(complex64(calc_eta_by_theta(k, abeli)) - etai) > tol):
        return - abeli - 0.5 * (1 + taufrom(k=k))
    else :
        return abeli


def calc_eta_by_theta(k, z):
    return 0.25 * complex(0, 1) * pi * ((jtheta(2, 0,  qfrom(k=k), 0)) ** 2) \
           * ((jtheta(4, 0,  qfrom(k=k), 0)) ** 2) \
           * (jtheta(3, 0,  qfrom(k=k), 0)) \
           * (jtheta(3, 2*z*pi,  qfrom(k=k), 0)) \
           / ( ((jtheta(1, z*pi,  qfrom(k=k), 0)) ** 2) * ((jtheta(3, z*pi,  qfrom(k=k), 0)) ** 2) )

def calc_mu(k, x1, x2, x3, zeta, abel):
    mu = []
    for i in range(0, 4, 1):
        mu.append(complex(
            0.25 *pi* ((jtheta(1, abel[i]*pi,  qfrom(k=k), 1) / (jtheta(1, abel[i]*pi,  qfrom(k=k), 0)) )
                       + (jtheta(3, abel[i]*pi,  qfrom(k=k), 1) / (jtheta(3, abel[i]*pi,  qfrom(k=k), 0)) )) \
            - x3 - (x2 + complex(0,1) *x1) * zeta[i]))
    return mu


def is_awc_multiple_root(k, x1, x2, x3):   # This will test if there are multiple roots; the analytic derivation assumes they are distinct
    K = complex64(ellipk(k**2))
    k1 = sqrt(1-k**2)
    tol = 0.001

    if ( (abs(4 * x1**2 * k1**2 - k**2 *( K**2 *k1**2 + 4* x2**2)) < tol) and x3==0):
        return True
    elif ( (abs(4 * x1**2 * k1**2 - K**2 *k1**2 + 4* x3**2 )< tol) and x2==0):
        return True

    return False

def is_awc_branch_point(k, x1, x2, x3):   # This will test if we get a branch point as a roots; these are numerically unstable
    K = complex64(ellipk(k**2))
    k1 = sqrt(1-k**2)
    tol = 0.001

    if ( (abs(k1 * x1- k*x2) < tol) and x3==0):
        return True


    return False


def higgs_squared(k, x1, x2, x3):

    if (is_awc_multiple_root(k, x1, x2, x3) ):
        return float(255)/float(256)

    if (is_awc_branch_point(k, x1, x2, x3) ):
        return float(255)/float(256)

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

    GNUM = grams(zeta, mu, [x1, x2, x3], k)

    inv_gram = matrix(GNUM).I

    higgs = phis(zeta, mu, [x1, x2, x3], k)

    # print GNUM, inv_gram, higgs # for testing

    return  -(trace(matmul( matmul(higgs, inv_gram),  matmul(higgs, inv_gram) )).real)/2


def higgs_squared_on_xy_plane(k_value, x0, x1, y0, y1, partition_size):

    x_step = (x1 - x0) / partition_size
    y_step = (y1 - y0) / partition_size


    points = []
    last = 0
    for k in range(0, partition_size):
        for i in range(0, partition_size):
            x = x0 + i * x_step
            y = y0 + k * y_step

            value = higgs_squared(k_value, x, y, 0.00)
            bucket_value = int(floor(256 * value))
            if(bucket_value > 255 or bucket_value < 0):
                print i, k, bucket_value
                bucket_value = 255
            points.append(bucket_value)

            last = bucket_value

    return points

def higgs_squared_on_xz_plane(k_value, x0, x1, z0, z1, partition_size):

    x_step = (x1 - x0) / partition_size
    z_step = (z1 - z0) / partition_size


    points = []
    last = 0
    for k in range(0, partition_size):
        for i in range(0, partition_size):
            x = x0 + i * x_step
            z = z0 + k * z_step

            value = higgs_squared(k_value, x, 0.00, z)
            bucket_value = int(floor(256 * value))
            if(bucket_value > 255 or bucket_value < 0):
                print i, k, bucket_value
                bucket_value = 255
            points.append(bucket_value)

            last = bucket_value

    return points



def write_point_to_file(points, filename):

    """

    :rtype : object
    """
    fo = open(os.path.expanduser("~/Desktop/numerical monopoles/testing_higgs/" + filename), 'wb')
    byteArray = bytearray(points)
    fo.write(byteArray)
    fo.close()


# t0 = time.time()
# print higgs_squared(0.8, 2.01, 0, 0)
# t1 = time.time()
# print str(t1-t0)

# print "%.8f"%  higgs_squared(0.8 , 1 , 0.5 , 0.04)



# fo = open(os.path.expanduser("~/Desktop/numerical monopoles/hwb_xyhiggs_k0"), 'w' )
# for i in range(0, 400, 1):
#     fo.write("%4.2f %15.9f\n"% ( (float(2*i) +2)/100, higgs_squared(0.8 ,  0.0, (float(2*i) +2)/100, 0.00 )    ))
# fo.close()


# for k in range(0, 5, 1):
#     for i in range(0, 5, 1):
#         x = 0.05 + i * 0.05
#         y = 0.05 + k * 0.05
#         print int(floor(256 * higgs_squared(0.001 , x, y, 0.0 ) ))
#
#
# print higgs_squared(0.8 , 6 , 7 , 0.0)
#
# print int(floor(256 * (1.0- higgs_squared(0.8 , 0.05 , 0.05 , 0.0))))
# #

# p = higgs_squared_on_xy_plane(0.0001, 0.05, 10.05, 0.05, 10.05,  200)   # k, x-initial x-final, y-initial, y-final, partition size=no points between initial final
#
# write_point_to_file(p , 'hwb_xyhiggs_k0001')


# p = higgs_squared_on_xz_plane(0.0001, 0.05, 10.05, 0.05, 10.05,  200)   # k, x-initial x-final, y-initial, y-final, partition size=no points between initial final
#
# write_point_to_file(p , 'hwb_xzhiggs_k0001')

# print higgs_squared(0.8, 0.4, 0.3, 0.001)

# k=0.13
# print k, sqrt( higgs_squared(k, complex64(ellipk(k**2))/2+0.005, 0.0, 0.0) )