import shapely.geometry as sg
import fiona
import uuid
import pdal
import json
import pprint
import numpy as np
import laspy
import time
import rtree
import scipy.spatial as sp
import startin
import json
import CGAL

# thinning function 
def thinning(array):
    thinning_factor = 2 #thinning factor
    new_array = array[::thinning_factor]
    return new_array

# getting AHN points and cellsize (cell/grid size is gonna be coded in a way that all parameters can be provided though Json/txt file)
def AHN_points():
    cellSize = 100  
    pc = laspy.file.File('C:/Users/khale/Google Drive/Geomatics_Temp/Geo1101/Interpolation/Data/LAS_amsterdam_out.las', mode='r')
    pts = np.vstack((pc.x, pc.y, pc.z)).transpose()
    return pts, cellSize 

# creation of the grid from the INPUT pointclound AND cellSize 
def creatGrid (points, cellSize):
    if type (points) == list:
        points = np.array(points)
    Min_XYZ = [int (np.amin(points[:, 0:1])-1), int (np.amin(points[:, 1:2])-1), int (np.amin(points[:, 2:])-1)]
    MAX_XYZ = [int (np.amax(points[:, 0:1]))  , int (np.amax(points[:, 1:2]))  , int (np.amax(points[:, 2:]))]
    
    grid_x = []
    grid_y = []
    a = Min_XYZ[0]
    while a < MAX_XYZ[0]:
        grid_x.append(a)
        a += cellSize

    b = Min_XYZ[1]
    while b < MAX_XYZ[1]:
        grid_y.append(b)
        b += cellSize

    grid = np.meshgrid(grid_x, grid_y)

    grid1 = grid[0].reshape((np.prod(grid[0].shape),))
    grid2 = grid[1].reshape((np.prod(grid[1].shape),))
    
    final_grid = []
    temp=[]
    grid_S = list(zip(grid1, grid2))
    for i in grid_S:
        temp.append((i[0], i[1],-9999))
        if i[0] == grid1[-1]:
            final_grid.append(temp)
            temp = []

    final_grid = np.asarray (final_grid)
    return final_grid, Min_XYZ[0], Min_XYZ[1]

# create and ascci file by providing the cellsize, the output file name/Directory x,y lowercorner  
def asci_creation(cellS, griPts, outPutas, xll, yll):
    nodata = -99999
    interp = open (outPutas, "w+")
    interp.write("NCOLS {0}\n".format(len(griPts[0])))
    interp.write("NROWS {0}\n".format(len(griPts)))
    interp.write("XLLCORNER {0}\n".format(xll))
    interp.write("YLLCORNER {0}\n".format(yll))
    interp.write("CELLSIZE {0}\n".format(cellS))
    interp.write("NODATA_VALUE {0}\n".format(nodata))

    #loop to write z to asci cells
    # count = 0
    for i in range (len (griPts)):
        for j in range(len(griPts[0])):
            print (griPts[i,j][2])
            interp.write("{0} ".format(griPts[i,j][2]))
        interp.write('\n')
    interp.close()

def main():
    ascc_output     = "test_output.asc"
    pts, cellSize   = AHN_points()
    grid, xll, yll  = creatGrid(pts, cellSize)
    asci_creation (cellSize, grid, ascc_output,xll, yll )
    # print (grid)

if __name__ == '__main__':
    main()

print ("Check")

