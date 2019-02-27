from numpy import sqrt, int, arange
import os
import numpy
import copy
from PIL import Image
import data

def reconstruct_2d(data):
    n = int(sqrt(data.size))
    points = []
    for i in range(0, n):
        points.append([])
    for i in range(0, data.size):
        y = i / n
        points[y].append(data[i])
    return points
    
def reflect_symmetries(positive_quadrant):
    bottom_right_quadrant = positive_quadrant
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

def png_from(pth):
    dat = unsmoothed_data(pth)
    flattened = [x for sub in dat for x in sub]

    im = Image.new('L', (60, 60), 'black') 
    im.putdata(flattened)

    im.save('./out.png')

def load_data():
    zrange = arange(0.025, 2.975, 0.05)
    krange = arange(0.01, 0.99, 0.01)       # This is Peter's default. Will try hardfix
#   krange = arange(0.95,0.96)

    im = Image.new('L', (60 * len(krange), 60 * len(zrange)), 'white')
    ims = []

    for k in krange:
        zres = []
        zim = Image.new( 'L', (60, 60 * len(zrange)), 'black')
        for z in zrange:
            pth = './python_smoothed/k=%.02f/xy_%.02f_%.03f' % (k, k, z)
            print pth
            dat = data.unsmoothed_data(pth)
            flattened = [x for sub in dat for x in sub]
            flattened = numpy.nan_to_num(flattened)
            # value is in range 0 - 1.5
            scaled = [x / 1.5 for x in flattened]
            scaled = [int(x * 255) for x in scaled]
            zres += scaled

        zim.putdata(zres)
        ims.append(zim)

    x_offset = 0
    for zslice in ims:
        im.paste(zslice, (x_offset, 0))
        x_offset += zslice.size[0]
       
    im.save('./out.png')

#png_from('./python_results/k=0.10/xy_0.1_0.025-3.025_0.025-3.025_0.525_60')
load_data()
