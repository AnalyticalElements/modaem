#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple


SurferGrid = namedtuple(
    "SurferGrid", ["nx", "ny", "xmin", "xmax", "ymin", "ymax", "grid"]
)


def read_grd(filename):
    with open(filename) as infile:
        fmt = str(infile.readline())
        if fmt != "DSAA\n":
            raise Exception("Wrong file format...")
        nx, ny = np.array(infile.readline().split(), dtype=np.int)
        xmin, xmax = np.array(infile.readline().split(), dtype=np.float)
        ymin, ymax = np.array(infile.readline().split(), dtype=np.float)
        grid_min, grid_max = np.array(infile.readline().split(), dtype=np.float)
    grid = np.loadtxt(filename, skiprows=5)

    return SurferGrid(nx, ny, xmin, xmax, ymin, ymax, grid)


def plot(filename, levels=30):
    surfer = read_grd(filename)
    plt.figure()
    plt.contour(surfer.grid, levels=levels)
    plt.show()

if __name__ == "__main__":
    plot(sys.argv[1])
