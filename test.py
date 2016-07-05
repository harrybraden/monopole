__author__ = 'hwb'

import time

from energy_density import calc_zeta, calc_eta, calc_abel, calc_mu
from python_expressions.dmus import dmus
from python_expressions.dzetas import dzetas
from python_expressions.ddmus import ddmus
from python_expressions.ddzetas import ddzetas
from python_expressions.ddphis111 import ddphis111

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

    print str(t1-t0)
    print str(t2-t1)
    print str(t3-t2)

    return A # ddphis111(zeta, mu, DM, DZ, DDM,  DDZ, [x1, x2, x3], k)

t4 = time.time()
B  = test_timing(.8, 1.5, 0.5, 0.2)
t5 = time.time()

print B
print str(t5-t4)