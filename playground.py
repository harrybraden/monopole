__author__ = 'hwb'
from numpy import roots, complex, complex64, mat, dot, trace, pi, sqrt, sum, trace, linalg, matmul, array, matrix, conj
from cmath import exp
import time
from mpmath import ellipk, ellipe, j, taufrom, jtheta, qfrom, ellipf, asin, mfrom
from numpy import roots, complex64, conj, pi, sqrt, sum, trace, linalg, matmul, array
import time
import math


from python_expressions.dexp import dexp
from python_expressions.dmus import dmus
from python_expressions.dzetas import dzetas
from python_expressions.ddmus import ddmus
from python_expressions.ddzetas import ddzetas

from python_expressions.grams import grams
from python_expressions.dgrams1 import dgrams1
from python_expressions.dgrams2 import dgrams2
from python_expressions.dgrams3 import dgrams3

from python_expressions.phis import phis
from python_expressions.dphis111 import dphis111
from python_expressions.dphis112 import dphis112
# from python_expressions.dphis121 import dphis121
from python_expressions.dphis122 import dphis122
from python_expressions.dphis211 import dphis211
from python_expressions.dphis212 import dphis212
from python_expressions.dphis221 import dphis221
from python_expressions.dphis222 import dphis222
from python_expressions.dphis311 import dphis311
from python_expressions.dphis312 import dphis312
# from python_expressions.dphis321 import dphis321
from python_expressions.dphis322 import dphis322
from python_expressions.ddgrams111 import ddgrams111
from python_expressions.ddgrams112 import ddgrams112
from python_expressions.ddgrams121 import ddgrams121
from python_expressions.ddgrams122 import ddgrams122
from python_expressions.ddgrams211 import ddgrams211
from python_expressions.ddgrams212 import ddgrams212
# from python_expressions.ddgrams221 import ddgrams221
from python_expressions.ddgrams222 import ddgrams222
from python_expressions.ddgrams311 import ddgrams311
from python_expressions.ddgrams312 import ddgrams312
# from python_expressions.ddgrams321 import ddgrams321
from python_expressions.ddgrams322 import ddgrams322
from python_expressions.ddphis111 import ddphis111
from python_expressions.ddphis112 import ddphis112
# from python_expressions.ddphis121 import ddphis121
from python_expressions.ddphis122 import ddphis122
from python_expressions.ddphis211 import ddphis211
from python_expressions.ddphis212 import ddphis212
# from python_expressions.ddphis221 import ddphis221
from python_expressions.ddphis222 import ddphis222
from python_expressions.ddphis311 import ddphis311
from python_expressions.ddphis312 import ddphis312
# from python_expressions.ddphis321 import ddphis321
from python_expressions.ddphis322 import ddphis322

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

    energy_density = -(ed1 + ed2 + ed3).real

    return energy_density






k = 0.2

x1 = 3.2
x2 = 1.3
x3 = 1.4

t1 = time.time()
A  = test(k , x1, x2, x3)
t2 = time.time()

print A

print str(t2-t1)

# print test(k , x1, x2, x3)


