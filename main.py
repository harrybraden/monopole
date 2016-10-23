__author__ = 'hwb'

from energy_density import energy_density, energy_density_on_xy_plane, write_point_to_file
from higgs_squared import higgs_squared
import argparse
import time

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('k', type=float)
    parser.add_argument('x_initial', type=float)
    parser.add_argument('z_final', type=float)
    parser.add_argument('z_initial', type=float)
    parser.add_argument('y_final', type=float)
    parser.add_argument('y', type=float)
    parser.add_argument('partition_size', type=int)
    args = parser.parse_args()

    print "Rendering XY[k=%s, x initial=%s, z final=%s, z initial=%s, y final=%s, y=%s, partition=%s]" % (args.k, args.x_initial, args.z_final, args.z_initial, args.y_final, args.y, args.partition_size)

    # k, x-initial z-final, z-initial, y-final, y, partition size=no points between initial final
    p = energy_density_on_xy_plane(args.k, args.x_initial, args.z_final, args.z_initial, args.y_final, args.y, args.partition_size)
    write_point_to_file(p , 'example_xy_%s_%s_%s' % (args.k, args.x_initial, args.y_final) )
