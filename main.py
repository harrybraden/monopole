__author__ = 'hwb'
from mpmath import ellipk, j, taufrom, jtheta, qfrom, ellipf, asin, mfrom
from numpy import roots, complex64, conj, pi, sqrt, sum
import mpmath
import time
import matrices
import os
import laplace



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


def calc_phi_squared(k, x1, x2, x3):
    t0= time.time()
    zeta = calc_zeta(k ,x1, x2, x3)
    t1 = time.time()
    # print "zeta: " + str(t1-t0)
    eta = calc_eta(k, x1, x2, x3)
    t2 = time.time()
    # print "eta: " + str(t2-t1)
    abel = calc_abel(k, zeta, eta)
    t3 = time.time()
    # print "abel: " + str(t3-t2)
    mu = calc_mu(k, x1, x2, x3, zeta, abel)
    t4 = time.time()
    # print "mu: " + str(t4-t3)

    result =  matrices.HIGGSTRACE(map(lambda z:complex(z), zeta), mu, [x1, x2, x3], k)
    t5 = time.time()
    # print "Higgs: " + str(t5-t4)
    # print "Total: " + str(t5-t0)
    return result.real

k = 0.8
x1 = 0.5
x2 = 0.0
x3 = 0.0

print "%.8f"%  calc_phi_squared(k ,x1, x2, x3).real
# print "%.8f"%  calc_phi_squared(k ,6, 6, 6).real


# phi = []
# for i in range(2, 600, 2):
#     phi.append(calc_phi_squared(0.8, float(i)/100, 0, 0).real)
# lap_phi = laplace(phi)*float( (50)**2 )
#
# fo = open(os.path.expanduser("~/Desktop/hwb_xlaplace"), 'w' )
# for i in range(0, len(lap_phi), 1):
#     fo.write("%4.2f %15.9f\n"% ( (float(2*i) +2)/100, lap_phi[i]))
# fo.close()

phi = []
for i in range(200, 600, 1):
    phi.append(calc_phi_squared(0.8, float(i)/100, 0, 0).real)
lap_phi = laplace(phi)*float( (100)**2 )

fo = open(os.path.expanduser("~/Desktop/hwb_xlaplaceR"), 'w' )
for i in range(0, len(lap_phi), 1):
    fo.write("%4.2f %15.9f\n"% ( (float(i)+200 )/100, lap_phi[i]))
fo.close()


print phi[100], phi[101], phi[102]
print lap_phi[100], lap_phi[101], lap_phi[102]


# t6 = time.time()
# fo = open(os.path.expanduser("~/Desktop/hwb_xdiag"), 'w' )
# for i in range(0,600,1):
#      fo.write("%4.2f %15.9f\n"% ((-5.99+i*.02)*sqrt(3), calc_phi_squared(k ,-5.99+i*.02, -5.99+i*.02, -5.99+i*.02).real) )
# fo.close()
# t7 = time.time()
#
# print t7-t6


zeta= calc_zeta(k, x1, x2, x3)
eta= calc_eta(k, x1, x2, x3)
abel= calc_abel(k, zeta, eta)
mu= calc_mu(k, x1, x2, x3, zeta, abel)

# print zeta
# print eta
# print abel
# print mu


# print sum(abel)
# print sum( abel+conj(abel) )
# print abel[0]+conj(abel)[2], abel[1]+conj(abel)[3]