__author__ = 'hwb'

from energy_density import energy_density, energy_density_on_xy_plane, write_point_to_file
from higgs_squared import higgs_squared
import argparse
import time

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('k', type=float)
    parser.add_argument('x0', type=float)
    parser.add_argument('x1', type=float)
    parser.add_argument('y0', type=float)
    parser.add_argument('y1', type=float)
    parser.add_argument('z', type=float)
    parser.add_argument('partition_size', type=int)
    args = parser.parse_args()

    print ("Rendering XY [k=%s, x0=%s, x1=%s, y0=%s, y1=%s, z=%s, partition=%s]" %
        (args.k, args.x0, args.x1, args.y0, args.y1, args.z, args.partition_size))

    #energy_density_on_xy_plane(k, x0, x1, y0, y1, z, partition_size):
    p = energy_density_on_xy_plane(args.k, args.x0, args.x1, args.y0, args.y1, args.z, args.partition_size)
    write_point_to_file(p , 'xy_%s_%s-%s_%s-%s_%s_%s' 
            % (args.k, args.x0, args.x1, args.y0, args.y1, args.z, args.partition_size) )
