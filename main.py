__author__ = 'hwb'
from energy_density import energy_density
from higgs_squared import higgs_squared
import time

t0 = time.time()
A  = energy_density(.8, 1.5, 0.0, 0.0)
t1 = time.time()

print A
print str(t1-t0)

t2 = time.time()
B  = higgs_squared(.8, 1.5, 0.0, 0.0)
t3 = time.time()

print B
print str(t3-t2)