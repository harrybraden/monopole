import os
import numpy
from numpy import sqrt
import copy
import array_tools
import smoothing_tools
from matplotlib.pyplot import * 


DIRECTORY = os.environ['MONOPOLE_OUTPUT'] or "~/Desktop/numerical monopoles/python_results/"
FILE = 'example_xy_2_3_30'
fo = open(os.path.expanduser(DIRECTORY + FILE), 'rb')
bytes = numpy.fromfile(fo, dtype=numpy.uint8)

shaped_list = array_tools.reconstruct_2d(bytes)
smoothed = smoothing_tools.smooth_2d(shaped_list)
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

imgplot = imshow(full)
savefig(DIRECTORY + FILE + '.png')
