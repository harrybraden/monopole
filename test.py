__author__ = 'hwb'

import time

from energy_density import calc_zeta, calc_eta, calc_abel, calc_mu
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

t15 = time.time()
C  = test_timing(.8, 1.5, 0.5, 0.2)
t16 = time.time()

print C
# print str(t16-t15)