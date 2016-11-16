# from __future__ import print_function

#  This file will calculate the energy_density for a point (x1, x2, x3) of space and a parameter k (between 0 and 1)
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
#  The energy density is one such function and we have the gram matrix (grams) and higgs field phis, both 2x2 matrices and their first and second derivatives dgrams, ddgrams,
#  dphis, ddphis
#
#
#
#
#

__author__ = 'hwb'

from numpy import roots, complex, complex64, complex128, mat, dot, trace, pi, sqrt, sum, trace, linalg, matmul, array, matrix, conj,  matmul, floor
from cmath import exp
import time
from mpmath import ellipk, ellipe, j, taufrom, jtheta, qfrom, ellipf, asin, mfrom
# from numpy import roots, complex64, conj, pi, sqrt, sum, trace, linalg, array
import time
import math
import os
import sys

maxint = 256                     # This determines the digits being returned
maxintr = maxint -1


from python_expressions.dexp import dexp
from python_expressions.dmus import dmus
from python_expressions.dzetas import dzetas
from python_expressions.ddmus import ddmus
from python_expressions.ddzetas import ddzetas

from python_expressions.grams import grams
from modified_expressions.dgrams1 import dgrams1
from modified_expressions.dgrams2 import dgrams2
from modified_expressions.dgrams3 import dgrams3

from python_expressions.phis import phis
from modified_expressions.dphis111 import dphis111
from modified_expressions.dphis112 import dphis112
# from modified_expressions.dphis121 import dphis121
from modified_expressions.dphis122 import dphis122
from modified_expressions.dphis211 import dphis211
from modified_expressions.dphis212 import dphis212
# from modified_expressions.dphis221 import dphis221
from modified_expressions.dphis222 import dphis222
from modified_expressions.dphis311 import dphis311
from modified_expressions.dphis312 import dphis312
# from modified_expressions.dphis321 import dphis321
from modified_expressions.dphis322 import dphis322
from modified_expressions.ddgrams111 import ddgrams111
from modified_expressions.ddgrams112 import ddgrams112
# from modified_expressions.ddgrams121 import ddgrams121
from modified_expressions.ddgrams122 import ddgrams122
from modified_expressions.ddgrams211 import ddgrams211
from modified_expressions.ddgrams212 import ddgrams212
# from modified_expressions.ddgrams221 import ddgrams221
from modified_expressions.ddgrams222 import ddgrams222
from modified_expressions.ddgrams311 import ddgrams311
from modified_expressions.ddgrams312 import ddgrams312
# from modified_expressions.ddgrams321 import ddgrams321
from modified_expressions.ddgrams322 import ddgrams322
from modified_expressions.ddphis111 import ddphis111
from modified_expressions.ddphis112 import ddphis112
# from modified_expressions.ddphis121 import ddphis121
from modified_expressions.ddphis122 import ddphis122
from modified_expressions.ddphis211 import ddphis211
from modified_expressions.ddphis212 import ddphis212
# from modified_expressions.ddphis221 import ddphis221
from modified_expressions.ddphis222 import ddphis222
from modified_expressions.ddphis311 import ddphis311
from modified_expressions.ddphis312 import ddphis312
# from modified_expressions.ddphis321 import ddphis321
from modified_expressions.ddphis322 import ddphis322

# def eprint(*args, **kwargs):
#         print(*args, file=sys.stderr, **kwargs)

def quartic_roots(k, x1, x2, x3):
    K = complex128(ellipk(k**2))

    e0 = complex128((x1*j - x2)**2 + .25 * K**2)
    e1 = complex128(4*(x1*j-x2)*x3)
    e2 = complex128(4*(x3**2) - 2 * (x1**2) - 2 * (x2**2) + (K**2) * (k**2 - 0.5))
    e3 = complex128(4*x3*(x2 + j*x1))
    e4 = complex128(x2**2 - x1**2 + 2*j*x1*x2 + 0.25*K**2)

    return roots([e4, e3, e2, e1, e0])


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
    tol = 0.001

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
    tol1 = 0.05
    tol2 = 0.05

    # Two smoothings here. This one is the better

    if ( (abs(4 * x1**2 * k1**2 - k**2 *( K**2 *k1**2 + 4* x2**2)) < tol1) and abs(x3)<tol2):
        return True
    elif ( (abs(4 * x1**2 * k1**2 - K**2 *k1**2 + 4* x3**2 )< tol1) and abs(x2)<tol2):
        return True

    # Second smoothing

    # tol = 0.01
    #
    # if ( (abs(4 * x1**2 * k1**2 - k**2 *( K**2 *k1**2 + 4* x2**2)) < tol) and x3==0):
    #     return True
    # elif ( (abs(4 * x1**2 * k1**2 - K**2 *k1**2 + 4* x3**2 )< tol) and x2==0):
    #     return True

    return False

def is_awc_branch_point(k, x1, x2, x3):   # This will test if we get a branch point as a roots; these are numerically unstable
    K = complex64(ellipk(k**2))
    k1 = sqrt(1-k**2)
    tol = 0.001

    if ( (abs(k1 * x1- k * x2) < tol) and x3==0):
        return True


    return False

def energy_density(k, x1, x2, x3):   # If there is a multiple root or branch point a value of maxint-1 will be returned; maxint is globally set.

    if (is_awc_multiple_root(k, x1, x2, x3) ):
        return float(maxintr)/float(maxint)

    if (is_awc_branch_point(k, x1, x2, x3) ):
        return float(maxintr)/float(maxint)

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

    inv_gram = matrix(GNUM).I

    higgs = phis(zeta, mu, [x1, x2, x3], k)

    DGS1 = dgrams1(zeta, mu, DM, DZ, x, k)

    DGS2 = dgrams2(zeta, mu, DM, DZ, x, k)

    DGS3 = dgrams3(zeta, mu, DM, DZ, x, k)


    # Using the hermiticity properties its faster to evaluate the matrix entries just once and so if we evaluate
    # the *12 element the *21 is minus the conjugate of this
    # DHS1 = mat([[ dphis111(zeta, mu, DM, DZ, [x1, x2, x3], k), dphis112(zeta, mu, DM, DZ, [x1, x2, x3], k)],
    #            [ dphis121(zeta, mu, DM, DZ, [x1, x2, x3], k), dphis122(zeta, mu, DM, DZ, [x1, x2, x3], k)]])

    DH112 = dphis112(zeta, mu, DM, DZ, [x1, x2, x3], k)

    DHS1 = mat([[ dphis111(zeta, mu, DM, DZ, [x1, x2, x3], k), DH112   ],
                [ -conj(DH112), dphis122(zeta, mu, DM, DZ, [x1, x2, x3], k)]])

    DH212 = dphis212(zeta, mu, DM, DZ, [x1, x2, x3], k)

    DHS2 = mat([[ dphis211(zeta, mu, DM, DZ, [x1, x2, x3], k), DH212   ],
                [ -conj(DH212), dphis222(zeta, mu, DM, DZ, [x1, x2, x3], k)]])

    DH312 = dphis312(zeta, mu, DM, DZ, [x1, x2, x3], k)

    DHS3 = mat([[ dphis311(zeta, mu, DM, DZ, [x1, x2, x3], k), DH312   ],
                [ -conj(DH312), dphis322(zeta, mu, DM, DZ, [x1, x2, x3], k)]])

    DDGS112 = ddgrams112(zeta, mu, DM, DZ, DDM, DDZ, [x1, x2, x3], k)

    DDGS1 = mat([[ ddgrams111(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k), DDGS112 ],
                 [ -conj(DDGS112), ddgrams122(zeta, mu, DM, DZ, DDM, DDZ, [x1, x2, x3], k)]])

    DDGS212 = ddgrams212(zeta, mu, DM, DZ, DDM, DDZ, [x1, x2, x3], k)

    DDGS2 = mat([[ ddgrams211(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k), DDGS212 ],
                 [ -conj(DDGS212) , ddgrams222(zeta, mu, DM, DZ, DDM, DDZ, [x1, x2, x3], k)]])

    DDGS312 = ddgrams312(zeta, mu, DM, DZ, DDM, DDZ, [x1, x2, x3], k)

    DDGS3 = mat([[ ddgrams311(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k), DDGS312 ],
                 [ -conj(DDGS312), ddgrams322(zeta, mu, DM, DZ, DDM, DDZ, [x1, x2, x3], k)]])

    DDHS111 = ddphis111(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)
    DDHS112 = ddphis112(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)
    # DDHS121 = ddphis121(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)
    DDHS122 = ddphis122(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

    DDHS1 = mat( [[DDHS111, DDHS112], [ -conj(DDHS112),DDHS122]])

    DDHS211 = ddphis211(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)
    DDHS212 = ddphis212(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)
    # DDHS221 = ddphis221(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)
    DDHS222 = ddphis222(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

    DDHS2 = mat( [[DDHS211, DDHS212], [-conj(DDHS212),DDHS222]])

    DDHS311 = ddphis311(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)
    DDHS312 = ddphis312(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)
    # DDHS321 = ddphis321(zeta, mu,DM, DZ, DDM,  DDZ, [x1, x2, x3], k)
    DDHS322 = ddphis322(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

    DDHS3 = mat( [[DDHS311, DDHS312], [-conj(DDHS312),DDHS322]])

    ed1 = trace(matmul( matmul(DDHS1, inv_gram) -2* matmul( matmul(DHS1 , inv_gram), matmul(DGS1, inv_gram)) \
                        + matmul(higgs, matmul( 2* matmul(matmul(inv_gram,DGS1), matmul(inv_gram, DGS1)), inv_gram)  - matmul(matmul(inv_gram, DDGS1), inv_gram)),
                        matmul(higgs,inv_gram)) ) \
          + trace( matmul(matmul(DHS1, inv_gram) - matmul(matmul(higgs, inv_gram), matmul(DGS1, inv_gram)),
                          matmul(DHS1, inv_gram) - matmul(matmul(higgs, inv_gram), matmul(DGS1, inv_gram))))

    ed2 = trace(matmul( matmul(DDHS2, inv_gram) -2* matmul( matmul(DHS2 , inv_gram), matmul(DGS2, inv_gram)) \
                        + matmul(higgs, matmul( 2* matmul(matmul(inv_gram,DGS2), matmul(inv_gram, DGS2)), inv_gram)  - matmul(matmul(inv_gram, DDGS2), inv_gram)),
                        matmul(higgs,inv_gram)) ) \
          + trace( matmul(matmul(DHS2, inv_gram) - matmul(matmul(higgs, inv_gram), matmul(DGS2, inv_gram)),
                          matmul(DHS2, inv_gram) - matmul(matmul(higgs, inv_gram), matmul(DGS2, inv_gram))))

    ed3 = trace(matmul( matmul(DDHS3, inv_gram) -2* matmul( matmul(DHS3 , inv_gram), matmul(DGS3, inv_gram)) \
                        + matmul(higgs, matmul( 2* matmul(matmul(inv_gram,DGS3), matmul(inv_gram, DGS3)), inv_gram)  - matmul(matmul(inv_gram, DDGS3), inv_gram)),
                        matmul(higgs,inv_gram)) ) \
          + trace( matmul(matmul(DHS3, inv_gram) - matmul(matmul(higgs, inv_gram), matmul(DGS3, inv_gram)),
                          matmul(DHS3, inv_gram) - matmul(matmul(higgs, inv_gram), matmul(DGS3, inv_gram))))

    # energy_density = -(ed1 + ed2 + ed3).real

    return  -(ed1 + ed2 + ed3).real

def energy_density_at_origin(k):
    K = complex64(ellipk(k**2))
    E = complex64(ellipe(k**2))
    k1 = sqrt(1-k**2)

    A = 32*(k**2 *(-K**2 * k**2 +E**2-4*E*K+3* K**2 + k**2)-2*(E-K)**2)**2/(k**8 * K**4 * k1**2)

    return A.real

def energy_density_on_xy_plane(k, x0, x1, y0, y1, z, partition_size):  # If this falls outside of [0,1) an value of maxint-1 will be returned; maxint is globally set.

    x_step = (x1 - x0) / partition_size
    y_step = (y1 - y0) / partition_size

    points = []

    for j in xrange(0, partition_size):
        # if j % 10 == 0 and j > 0:
        #     eprint("- rendered %s lines..." % j)

        for i in xrange(0, partition_size):
            x = x0 + i * x_step
            y = y0 + j * y_step

            value = energy_density(k, x, y, z)
            bucket_value = int(floor(maxint * value))
            if(bucket_value > maxintr or bucket_value < 0):
                print i, j, bucket_value
                bucket_value = maxintr
            points.append(bucket_value)



    return points

def energy_density_on_yz_plane(k, y0, y1, z0, z1, x, partition_size):   # If this falls outside of [0,1) an value of maxint-1 will be returned; maxint is globally set.

    y_step = (y1 - y0) / partition_size
    z_step = (z1 - z0) / partition_size

    points = []
    for j in range(0, partition_size):
        for i in range(0, partition_size):
            y = y0 + i * y_step
            z = z0 + j * z_step

            value = energy_density(k, x, y, z)
            bucket_value = int(floor(maxint * value))
            if(bucket_value > maxintr or bucket_value < 0):
                print i, j, bucket_value
                bucket_value = maxintr
            points.append(bucket_value)


    return points

def energy_density_on_xz_plane(k, x0, x1, z0, z1, y, partition_size):   # If this falls outside of [0,1) an value of maxint-1 will be returned; maxint is globally set.

    x_step = (x1 - x0) / partition_size
    z_step = (z1 - z0) / partition_size


    points = []
    for j in range(0, partition_size):
        for i in range(0, partition_size):
            x = x0 + i * x_step
            z = z0 + j * z_step

            value = energy_density(k, x, y, z)
            bucket_value = int(floor(maxint * value))
            if(bucket_value > maxintr or bucket_value < 0):
                print i, j, bucket_value
                bucket_value = maxintr
            points.append(bucket_value)


    return points



def test_timing(k, x1, x2, x3):

    t0 = time.time()
    zeta = calc_zeta(k ,x1, x2, x3)
    eta = calc_eta(k, x1, x2, x3)
    abel = calc_abel(k, zeta, eta)
    mu = calc_mu(k, x1, x2, x3, zeta, abel)
    x=[x1,x2,x3]


    t1 = time.time()

    K = complex64(ellipk(k**2))
    xp = x[0]+complex(0,1)*x[1]
    xm = x[0]-complex(0,1)*x[1]


    DM = dmus(zeta, x, k)
    DZ = dzetas(zeta, x,k)
    DDM = ddmus(zeta, x, k)
    DDZ = ddzetas(zeta, x,k)

    t2 =  time.time()

    A = ddphis111(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

    t3= time.time()

    B = ddphis222(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

    t4= time.time()

    C = ddgrams211(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

    t5= time.time()
    return A, B, C

def write_point_to_file(points, filename):
    """
    :rtype : object
    """
    fo = open(os.path.expanduser(filename), 'wb')
    byteArray = bytearray(points)
    fo.write(byteArray)
    fo.close()



# t15 = time.time()
# D  = test_timing(.8, 1.5, 0.5, 0.2)
# t16 = time.time()
#
# print D
# print str(t16-t15)
#
# t4 = time.time()
# A  = energy_density(.8, 1.5, 0.7, 0.3)
# t5 = time.time()
#
# print A
# print str(t5-t4)


# print order_roots(quartic_roots(0.8, 1.0, 0, 2.35))


# print quartic_roots(0.8, 1 , 0, 2.35)
