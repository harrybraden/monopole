from numpy import sqrt, int, arange
import numpy as np
import os
import numpy
import copy
import smoothing_tools
import mcubes
import argparse
import scipy.ndimage
import data

def print_obj(points, faces):
    print "# OBJ file"
    for v in points:
        print ("v %d %d %d" % (v[0], v[1], v[2]))

    for f in faces:
        print ("f %d %d %d" % (f[0] + 1, f[1] + 1, f[2] + 1))

def create_mesh(k, intensity):
    volume = data.load_data(k)
    vertices, triangles = mcubes.marching_cubes(volume, intensity)
    print_obj(vertices, triangles)
    # TODO use Catmull-Clarke to smooth voxel mesh
    print "#", k, intensity



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('k', type=float)
    parser.add_argument('threshold', type=float)
    args = parser.parse_args()
    create_mesh(args.k, args.threshold)
