import os
import numpy
from numpy import sqrt
import copy
import array_tools
import smoothing_tools
import matplotlib.pyplot as plt
import argparse

def load(path):
    fo = open(os.path.expanduser(path), 'rb')
    bytes = numpy.fromfile(fo, dtype=numpy.uint8)
    shaped_list = array_tools.reconstruct_2d(bytes)
    return shaped_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    args = parser.parse_args()
    print "Rendering %s" % args.file

    data = load(args.file)
    #smoothed = smoothing_tools.smooth_2d(data)
    reflected = array_tools.reflect_symmetries(data)

    plt.title(args.file)
    imgplot = plt.imshow(reflected, cmap=plt.get_cmap('viridis'), vmin=0, vmax=255)
    plt.savefig(args.file + '.png')
