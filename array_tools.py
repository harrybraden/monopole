from numpy import *
import os
from os import listdir
from os.path import isfile, join
import numpy
import copy
import smoothing_tools

def flatten_2d(data):
    return None


def flatten_3d(data):
    return None


def reconstruct_2d(data):
    n = int(sqrt(data.size))
    points = []
    for i in range(0, n):
        points.append([])
    for i in range(0, data.size):
        y = i / n
        points[y].append(data[i])
    return points

def reconstruct_3d(data):
    return None

def smoothed_image(data):
    file = "~/Documents/MonopoleMovie/python_converted/%s"  % data
    # file = "~/Desktop/numerical monopoles/exceptional/%s"  % data

    fo = open(os.path.expanduser(file), 'rb')

    bytes = numpy.fromfile(fo, dtype=numpy.uint8)

    shaped_list = reconstruct_2d(bytes)

    smoothed = smoothing_tools.smooth_2d(shaped_list)

    return reflect_symmetries(smoothed)

    
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


def unsmoothed_image(data):
    file = "~/Desktop/numerical monopoles/python_results/%s"  % data
    # file = "~/Desktop/numerical monopoles/exceptional/%s"  % data

    fo = open(os.path.expanduser(file), 'rb')

    bytes = numpy.fromfile(fo, dtype=numpy.uint8)

    shaped_list = reconstruct_2d(bytes)

    return reflect_symmetries(shaped_list)


def load_slice(file):
    fo = open(file, 'rb')
    bytes = numpy.fromfile(fo, dtype=numpy.uint8)
    shaped_list = reconstruct_2d(bytes)
    return shaped_list

def get_z_value(file):
    parts = str.split(file, '_')
    return float(parts[len(parts) - 1])

def write_point_to_file(points, filename):
    """
    :rtype : object
    """
    fo = open(os.path.expanduser(filename), 'wb')
    byteArray = bytearray(points)
    fo.write(byteArray)
    fo.close()

def get_point_tuples(kk, e_min, e_max):
    k = "%1.2f" % kk
    points=[]

    x0 = 0.025
    xn = 3.025
    y0 = 0.025
    yn = 3.025
    n = 60
    dx = (xn - x0) / n
    dy = (yn - y0) / n

    directory = '/Users/hwb/Desktop/numerical monopoles/testing_energydensity'

    k_directory = directory + '/k=' + str(k) + '/'
    files = [f for f in listdir(k_directory) if isfile(join(k_directory, f))]
    for f in files:

        shaped_list = load_slice(k_directory + f)

        z = get_z_value(f)

        for j, xl in enumerate(shaped_list):
            for i, e in enumerate(xl):
                if (e >=e_min and  e <=e_max):
                    points.append(( x0 + i*dx, y0 + j*dy, z))
    return points


def get_max_value(k):

    k_directory = directory + '/k=' + k +'/'

    files = [f for f in listdir(k_directory) if isfile(join(k_directory, f))]

    max = 0
    for f in files:
        shaped_list = load_slice(k_directory + f)
        z = get_z_value(f)
        for i, xl in enumerate(shaped_list):
            for j, e in enumerate(xl):
                if (e > max and e != 255):
                    max = e
    return max

def get_exceptional_tuples(kk):
    k = "%1.2f" % kk
    exceptional = []

    x0 = 0.025
    xn = 3.025
    y0 = 0.025
    yn = 3.025
    n = 60
    dx = (xn - x0) / n
    dy = (yn - y0) / n

    directory = '/Users/hwb/Desktop/numerical monopoles/testing_energydensity'

    k_directory = directory + '/k=' + k +'/'
    files = [f for f in listdir(k_directory) if isfile(join(k_directory, f))]
    for f in files:

        shaped_list = load_slice(k_directory + f)

        z = get_z_value(f)

        for j, xl in enumerate(shaped_list):
            for i, e in enumerate(xl):
                if (e == 255):
                # if (e>230):
                    exceptional.append((x0 + i*dx, y0 + j*dy, z))
    return exceptional

# print get_exceptional_tuples(0.80)
# print get_point_tuples(0.80,253,256)

def add_symmetric_points(original_points):
    points = []
    for t in original_points:
        points.append((t[0], t[1], t[2]))
        points.append((-t[0], t[1], t[2]))
        points.append((-t[0], -t[1], t[2]))
        points.append((t[0], -t[1], t[2]))
        points.append((t[0], t[1], -t[2]))
        points.append((-t[0], t[1], -t[2]))
        points.append((-t[0], -t[1], -t[2]))
        points.append((t[0], -t[1], -t[2]))
    return points


def plot_energy_density(k, lower_range, upper_range):
    points = add_symmetric_points(get_point_tuples(k, lower_range, upper_range))
    # points = add_symmetric_points(get_exceptional_tuples())

    x = map(lambda t: t[0], points)
    y = map(lambda t: t[1], points)
    z = map(lambda t: t[2], points)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    ax.set_xlim3d(-2.0,2.0)
    ax.set_ylim3d(-2.0,2.0)
    ax.set_zlim3d(-2.0,2.0)

    ax.view_init(elev=18., azim=10)

    ax.scatter(x, y, z, s=1)


def write_point_tuples(k, e_min, e_max):
    points=[]

    x0 = 0.025
    xn = 3.025
    y0 = 0.025
    yn = 3.025
    n = 60
    dx = (xn - x0) / n
    dy = (yn - y0) / n

    k_directory = directory + '/k=' + k + '/'
    files = [f for f in listdir(k_directory) if isfile(join(k_directory, f))]
    for f in files:

        shaped_list = load_slice(k_directory + f)

        z = get_z_value(f)

        for i, xl in enumerate(shaped_list):
            for j, e in enumerate(xl):
                if (e > e_min and  e < e_max):
                    points.append((x0 + i*dx, y0 + j*dy, z))
    return points

def plot_energy_density(k, lower_range, upper_range):
    points = add_symmetric_points(get_point_tuples(k, lower_range, upper_range))
    # points = add_symmetric_points(get_exceptional_tuples())

    x = map(lambda t: t[0], points)
    y = map(lambda t: t[1], points)
    z = map(lambda t: t[2], points)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    ax.set_xlim3d(-2.0,2.0)
    ax.set_ylim3d(-2.0,2.0)
    ax.set_zlim3d(-2.0,2.0)

    ax.view_init(elev=18., azim=10)

    ax.scatter(x, y, z, s=1)


def plot_exceptional(k):
    points = add_symmetric_points(get_exceptional_tuples(k))

    x = map(lambda t: t[0], points)
    y = map(lambda t: t[1], points)
    z = map(lambda t: t[2], points)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    ax.set_xlim3d(-0.5,0.5)
    ax.set_ylim3d(-0.8,0.8)
    ax.set_zlim3d(-2.0,2.0)

    ax.view_init(elev=18., azim=90)

    ax.scatter(x, y, z, s=1)


def write_point_tuples(k, e_min, e_max):
    points = add_symmetric_points(get_point_tuples(k, e_min, e_max))


    write_file = '/Users/hwb/Desktop/numerical monopoles/formaple/all_' + str(k)  +'.txt'

    # {[1, 2, 3],[3,5,6],}...

    fo = open(os.path.expanduser(write_file), 'w+')
    # fo.write('{')
    for i, point in enumerate(points):
        # point_string = '[' + str(point[0]) +','+ str(point[1]) + ',' + str(point[2]) + '] '
        point_string =  str(point[0]) +' '+ str(point[1]) + ' ' + str(point[2]) +'\n'
        if (i != (len(points) - 1)):
            point_string = point_string  #+ ','
        # print point_string
        fo.write(point_string)

    # fo.write('}')
    fo.close()

# write_point_tuples("0.95",60,80)



