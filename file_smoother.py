from array import array
import os
from numpy import mean
from math import floor

from os import listdir
from os.path import isfile, join

def smooth_value(i, data):
    good_points = []
    j = floor(i/60)
    if (data[i-2] > 0 and j == floor((i-2)/60)):
        good_points.append(data[i-2])
    if (data[i-1] > 0 and j == floor((i-1)/60)):
        good_points.append(data[i-1])
    if (data[i+1] > 0 and j == floor((i+1)/60)):
        good_points.append(data[i+1])
    if (data[i+2] > 0 and j == floor((i+2)/60)):
        good_points.append(data[i+2])

    return mean(good_points)

def convert_file(filename, converted_filename):
    input_file = open(os.path.expanduser(filename), 'r')
    data = array('d')
    data.read(input_file, 3600)
    data.tolist()
    smoothed_data = []

    for i, d in enumerate(data):

        smoothed_value = d
        if (d == 0 or i % 60 == 0):
            smoothed_value = smooth_value(i, data)
            # print 'smoothed ' + str(d) + ' to ' + str(smoothed_value)

        smoothed_data.append(smoothed_value)


    fo = open(os.path.expanduser(converted_filename), 'wb')
    fo.write(array('d', smoothed_data))
    fo.close()



directory = os.path.expanduser('~/Documents/MonopoleMovie/python_results')
outputdirectory= os.path.expanduser('~/Documents/MonopoleMovie/python_smoothed')




k_directories = [k for k in listdir(directory) if not isfile(k)]
for  k_directory in k_directories:
    files = [f for f in listdir(join(directory, k_directory)) if isfile(join(directory, k_directory, f))]
    for f in files:
        convert_file(join(directory, k_directory, f), join(outputdirectory, k_directory, f))
    print 'smoothed ' + k_directory

