from array import array
import os

from os import listdir
from os.path import isfile, join


def convert_file(filename, converted_filename):
    input_file = open(os.path.expanduser(filename), 'r')
    data = array('d')
    data.read(input_file, 3600)
    data.tolist()
    bucketed_data = []

    for d in data:
        bucketed_value = int(round(128.0 * d))
        if (bucketed_value >=256):
            bucketed_value = 255
        bucketed_data.append(bucketed_value)


    fo = open(os.path.expanduser(converted_filename), 'wb')
    byteArray = bytearray(bucketed_data)
    fo.write(byteArray)
    fo.close()



directory = os.path.expanduser('~/Documents/MonopoleMovie/python_results')
outputdirectory= os.path.expanduser('~/Documents/MonopoleMovie/python_converted')




k_directories = [k for k in listdir(directory) if not isfile(k)]
for  k_directory in k_directories:
    files = [f for f in listdir(join(directory, k_directory)) if isfile(join(directory, k_directory, f))]
    for f in files:
        convert_file(join(directory, k_directory, f), join(outputdirectory, k_directory, f))
    print 'converted ' + k_directory

