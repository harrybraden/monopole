__author__ = 'hwb'
from energy_density import energy_density
from higgs_squared import higgs_squared
import time

t0 = time.time()
A  = energy_density(.8, 1.5, 0.9, 0.3)
t1 = time.time()

print A
print str(t1-t0)

# t2 = time.time()
# B  = higgs_squared(.8, 1.5, 0.0, 0.0)
# t3 = time.time()
#
# print B
# print str(t3-t2)


# t0 = time.time()
# print higgs_squared(0.8, 2.01, 0, 0)
# t1 = time.time()
# print str(t1-t0)

# print "%.8f"%  higgs_squared(0.8 , 1 , 0.5 , 0.04)



# fo = open(os.path.expanduser("~/Desktop/numerical monopoles/hwb_testhiggs_2"), 'w' )
# for i in range(0, 400, 1):
#     fo.write("%4.2f %15.9f\n"% ( (float(2*i) +2)/100, higgs_squared(0.8 ,  0.0, (float(2*i) +2)/100, 0.00 )    ))
# fo.close()