from numpy import sqrt, int, arange
import os
import numpy
import copy
import smoothing_tools
from PIL import Image

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

def unsmoothed_data(pth):
    fo = open(os.path.expanduser(pth), 'rb')
    byte_data = numpy.fromfile(fo, dtype=numpy.uint8)
    return reconstruct_2d(byte_data)

def smoothed_data(pth):
    return smoothing_tools.smooth_2d(unsmoothed_data(pth))

def reflected_data(pth):
    return reflect_symmetries(smoothed_data(pth))

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
            # pth = './python_results/k=%.02f/xy_%s_0.025-3.025_0.025-3.025_%s_60' % (k, k, z)
            pth = './python_converted_scaled/k=%.02f/xy_%.02f_%s' % (k, k, z)
            print  pth
            dat = unsmoothed_data(pth)
            flattened = [x for sub in dat for x in sub]
            zres += flattened

        zim.putdata(zres)
        ims.append(zim)

    x_offset = 0
    for zslice in ims:
        im.paste(zslice, (x_offset, 0))
        x_offset += zslice.size[0]
       
    im.save('./out.png')

#png_from('./python_results/k=0.10/xy_0.1_0.025-3.025_0.025-3.025_0.525_60')
load_data()
