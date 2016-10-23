


def smooth_2d(grid):

    for i in range(0, len(grid)):
        for j in range(0, len(grid[i])):

            if (grid[i][j] == 255):
                neighbours = []
                neighbours.append(int(try_get(try_get(grid, i-1), j-1)))
                neighbours.append(int(try_get(try_get(grid, i-1), j)))
                neighbours.append(int(try_get(try_get(grid, i-1), j+1)))
                neighbours.append(int(try_get(try_get(grid, i), j-1)))
                neighbours.append(int(try_get(try_get(grid, i), j+1)))
                neighbours.append(int(try_get(try_get(grid, i+1), j-1)))
                neighbours.append(int(try_get(try_get(grid, i+1), j)))
                neighbours.append(int(try_get(try_get(grid, i+1), j+1)))

                neighbours = filter(lambda x: x!=255, neighbours)
                grid[i][j] = sum(neighbours) / len(neighbours)
    return grid

def try_get(array, index):
    try:
        return array[index]
    except:
        return 255

def smooth_3d(data):
    return None