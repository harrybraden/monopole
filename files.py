from array import array
import os


def write_floats(filename, floats):

    fo = open(os.path.expanduser(filename), 'wb')
    data = array('d', floats)
    data.write(fo)
    fo.close()

def read_floats(filename, size):
    input_file = open(os.path.expanduser(filename), 'r')
    data = array('d')
    data.read(input_file, size)   # 3600 here is the number of points in the file
    return data.tolist()


# floats = [2.26112841797, 2.25370947266]
#
# write_floats('tests', floats)
# data = read_floats('tests')
#
# for f in data:
#     print f
