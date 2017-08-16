from __future__ import print_function
from energy_density import energy_density, energy_density_on_xy_plane, write_point_to_file
import argparse
import time
import os
import sys

__author__ = 'hwb'


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('z', type=int)
    args = parser.parse_args()


    eprint ('true')

