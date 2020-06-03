### LAS READING AND RASTER PREPARATION ###

import math
import numpy as np
from laspy.file import File

def las_prepare(size, fpath, gnd_only = False):
    """Takes the filepath to an input LAS file and the
    desired output raster cell size. Reads the LAS file and outputs
    the ground points as a numpy array. Also establishes some
    basic raster parameters:
        - the extents
        - the resolution in coordinates
        - the coordinate location of the relative origin (bottom left)
    If called with gnd_only = True, it will ignore non-ground points,
    but this should optimally be done in the PDAL pipeline, not here.
    """
    in_file = File(fpath, mode = "r")
    header = in_file.header
    if gnd_only == True:
        in_np = np.vstack((in_file.raw_classification,
                           in_file.x, in_file.y, in_file.z)).transpose()
        in_np = in_np[in_np[:,0] == 2].copy()[:,1:]
    else: in_np = np.vstack((in_file.x, in_file.y, in_file.z)).transpose()
    extents = [[header.min[0], header.max[0]],
               [header.min[1], header.max[1]]]
    res = [math.ceil((extents[0][1] - extents[0][0]) / size),
           math.ceil((extents[1][1] - extents[1][0]) / size)]
    origin = [np.mean(extents[0]) - (size / 2) * res[0],
              np.mean(extents[1]) - (size / 2) * res[1]]
    return in_np, res, origin