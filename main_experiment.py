from __future__ import print_function
from energy_density import energy_density, energy_density_on_xy_plane, write_point_to_file
import argparse
import time
import os
import sys

__author__ = 'peterbraden'

# MONOPOLE_OUTPUT="~/Desktop/numerical monopoles/python_results/""

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('k', type=float)
    parser.add_argument('z', type=float)
    args = parser.parse_args()

    zval =-0.025 +0.05*args.z

    kt= args.k/100.0

    kk ="%1.2f" %kt
    
    eprint( kk)

    arglist = (eval(kk), 0.025, 3.025, 0.025, 3.025, zval, 60)
     
    eprint( arglist)

    eprint ("Rendering XY [k=%s, x0=%s, x1=%s, y0=%s, y1=%s, z=%s, partition=%s]" % arglist) 

    t_start = time.time()
#    p = energy_density_on_xy_plane(*arglist)
    t_end = time.time()

    # directory = os.environ['MONOPOLE_OUTPUT'] or './data'
    directory = "~/Desktop/numerical monopoles/python_results/"

    filename= (directory + 'xy_%s_%s-%s_%s-%s_%s_%s' % arglist)
#    write_point_to_file(p , filename)

    eprint ('-> finished [%ss], written to %s' % (str(t_end-t_start), filename))
