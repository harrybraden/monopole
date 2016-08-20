
from numpy import *

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