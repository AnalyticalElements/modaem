import numpy as np 

def read_grd(filename):
    with open(filename) as infile:
        fmt = infile.readline()
        nx, ny = np.array(infile.readline().split(), dtype=np.int)
        xmin, xmax = np.array(infile.readline().split(), dtype=np.float)
        ymin, ymax = np.array(infile.readline().split(), dtype=np.float)
        grid_min, grid_max = np.array(infile.readline().split(), dtype=np.float)
    grid = np.loadtxt(filename, skiprows=5)
         
    return nx, ny, xmin, xmax, ymin, ymax, grid

