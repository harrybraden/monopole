from __future__ import print_function
from energy_density import energy_density, energy_density_on_xy_plane, write_point_to_file
import argparse
import time
import os
import sys

__author__ = 'peterbraden'

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('k', type=float)
    parser.add_argument('xmax', type=float)
    parser.add_argument('ymax', type=float)
    parser.add_argument('z', type=float)
    parser.add_argument('partition_size', type=int)
    args = parser.parse_args()

    # Calculate x0, y0
    # Because we need to reflect across the axes, the minimum point is half
    # the distance between the points on the lattice.
    mp = 1./ (args.partition_size * 2)
    x0 = mp * args.xmax
    y0 = mp * args.ymax

    arglist = (args.k, x0, args.xmax, y0, args.ymax, args.z, args.partition_size)

    eprint ("Rendering XY [k=%s, x0=%s, x1=%s, y0=%s, y1=%s, z=%s, partition=%s]" % arglist) 

    t_start = time.time()
    p = energy_density_on_xy_plane(*arglist)
    t_end = time.time()

    directory = os.environ['MONOPOLE_OUTPUT'] or './data'
    filename= (directory + 'xy_%s_%s-%s_%s-%s_%s_%s' % arglist)
    write_point_to_file(p , filename)

    eprint ('-> finished [%ss], written to %s' % (str(t_end-t_start), filename))
