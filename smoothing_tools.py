


def smooth_2d(grid):

    for i in range(1, len(grid) - 1):
        for j in range(1, len(grid[i]) - 1):
            if (grid[i][j] == 255):
                grid[i][j] = (int(grid[i][j+1]) + int(grid[i][j-1]) + int(grid[i+1][j]) + int(grid[i-1][j]) )/ 4
    return grid


def smooth_3d(data):
    return None