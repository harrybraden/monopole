__author__ = 'hwb'
from numpy import roots, complex, complex64, mat, dot, trace, pi, sqrt, sum, trace, linalg, matmul, array, matrix
from cmath import exp
import time
from mpmath import ellipk, ellipe, j, taufrom, jtheta, qfrom, ellipf, asin, mfrom
from numpy import roots, complex64, conj, pi, sqrt, sum, trace, linalg, matmul, array
import time
import math
from generated_matrices.egram import egram
from generated_matrices.ephi import ephi
from generated_matrices.dgram1 import dgram1
from generated_matrices.dgram2 import dgram2
from generated_matrices.dgram3 import dgram3
from generated_matrices.ddgram1 import ddgram1
from generated_matrices.ddgram2 import ddgram2
from generated_matrices.ddgram3 import ddgram3
from generated_matrices.dphi1 import dphi1
from generated_matrices.dphi2 import dphi2
from generated_matrices.dphi3 import dphi3

from python_expressions.dexp import dexp
from python_expressions.dmus import dmus
from python_expressions.dzetas import dzetas
from python_expressions.grams import grams
from python_expressions.dgrams1 import dgrams1
from python_expressions.dgrams2 import dgrams2
from python_expressions.dgrams3 import dgrams3
from python_expressions.phis import phis
from python_expressions.dphis111 import dphis111
from python_expressions.dphis112 import dphis112
from python_expressions.dphis121 import dphis121
from python_expressions.dphis122 import dphis122
from python_expressions.dphis211 import dphis211
from python_expressions.dphis212 import dphis212
from python_expressions.dphis221 import dphis221
from python_expressions.dphis222 import dphis222
from python_expressions.dphis311 import dphis311
from python_expressions.dphis312 import dphis312
from python_expressions.dphis321 import dphis321
from python_expressions.dphis322 import dphis322

def quartic_roots(k, x1, x2, x3):
    K = complex64(ellipk(k**2))

    e0 = complex64((x1*j - x2)**2 + .25 * K**2)
    e1 = complex64(4*(x1*j-x2)*x3)
    e2 = complex64(4*(x3**2) - 2 * (x1**2) - 2 * (x2**2) + (K**2) * (k**2 - 0.5))
    e3 = complex64(4*x3*(x2 + j*x1))
    e4 = complex64(x2**2 - x1**2 + 2*j*x1*x2 + 0.25*K**2)

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



def test(k, x1, x2, x3):

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

    t0 = time.time()

    GNUM = egram(zeta, mu, [x1, x2, x3], k)

    t1 = time.time()

    # graminv = matrix(GNUM).I

    t2 = time.time()

    GNUM1 = grams(zeta, mu, [x1, x2, x3], k)

    t3 = time.time()

    t4 = time.time()
    higgs = ephi(zeta, mu, [x1, x2, x3], k)

    t5 = time.time()

    higgs1 = phis(zeta, mu, [x1, x2, x3], k)

    t6 = time.time()


    DH3 = dphi3(zeta, mu, [x1, x2, x3], k)

    t7 = time.time()

    DHS1 = mat([[ dphis111(zeta, mu, DM, DZ, [x1, x2, x3], k), dphis112(zeta, mu, DM, DZ, [x1, x2, x3], k)],
             [ dphis121(zeta, mu, DM, DZ, [x1, x2, x3], k), dphis122(zeta, mu, DM, DZ, [x1, x2, x3], k)]])

    DHS2 = mat([[ dphis211(zeta, mu, DM, DZ, [x1, x2, x3], k), dphis212(zeta, mu, DM, DZ, [x1, x2, x3], k)],
             [ dphis221(zeta, mu, DM, DZ, [x1, x2, x3], k), dphis222(zeta, mu, DM, DZ, [x1, x2, x3], k)]])

    t8 = time.time()

    DHS3 = mat([[ dphis311(zeta, mu, DM, DZ, [x1, x2, x3], k), dphis312(zeta, mu, DM, DZ, [x1, x2, x3], k)],
                [ dphis321(zeta, mu, DM, DZ, [x1, x2, x3], k), dphis322(zeta, mu, DM, DZ, [x1, x2, x3], k)]])

    t9 = time.time()

    # print str(t1-t0)
    # print str(t3-t2)
    # print str(t5-t4)
    # print str(t6-t5)
    print str(t7-t6)
    print str(t9-t8)

    return DH3, DHS3



    return DZ


k = 0.8

x1 = 0.6
x2 = 1.3
x3 = 2.0

print test(k , x1, x2, x3)


