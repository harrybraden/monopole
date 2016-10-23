import os
import numpy
from numpy import sqrt
import copy
import array_tools
import smoothing_tools
from matplotlib.pyplot import * 
import argparse

def load(path):
    fo = open(os.path.expanduser(path), 'rb')
    bytes = numpy.fromfile(fo, dtype=numpy.uint8)
    shaped_list = array_tools.reconstruct_2d(bytes)
    return shaped_list

def reflect(tl_quadrant):
    smoothed = smoothing_tools.smooth_2d(tl_quadrant)
    bottom_right_quadrant = smoothed
    top_right_quadrant = copy.deepcopy(bottom_right_quadrant)
    top_right_quadrant.reverse()
    right_half = top_right_quadrant + bottom_right_quadrant
    left_half = copy.deepcopy(right_half)

    for a in left_half:
        a.reverse()
        
    full = []
    for i in range(0, len(right_half)):
        full.append(None)
        full[i] = left_half[i] + right_half[i]

    return full


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    args = parser.parse_args()
    print "Rendering %s" % args.file

    data = load(args.file)
    reflected = reflect(data)

    imgplot = imshow(reflected, cmap=get_cmap('viridis'))
    savefig(args.file + '.png')
