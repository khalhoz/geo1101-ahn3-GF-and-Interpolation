import rasterio
import laspy 
import numpy as np

def read_file (directory):
	# function that readds raster (tiff) file and returns raster array and cellsize
    src = rasterio.open(directory)
    width = src.width
    height= src.height 
    bound = src.bounds
    raster_array = src.read(1)
    cell_size = ((src.transform * (src.width, src.height)) [0] - (src.transform * (0, 0))[0])/width, ((src.transform * (src.width, src.height)) [1] - (src.transform * (0, 0))[1])/height
    print (width, height, bound, len (raster_array), cell_size)
    return raster_array, cell_size

def read_PC_Data (directory): 
	# function that reads las file of PC_Data returns an array of the points and its classes
	# I am not sure the classes are actually needed here, you can see  
    pc = laspy.file.File(directory, mode='r')
    pc.classification
    pts = np.vstack((pc.x, pc.y, pc.z, pc.classification)).transpose()
    return pts

def cel_b_box(cell, sellsize): 
	# function that return bounding box of the cell (the cell itself)
	lowl = cell[0] - sellsize[0]/2, cell[1] - sellsize[1]/2
	lowr = cell[0] + sellsize[0]/2, cell[1] - sellsize[1]/2
	Up_r = cell[0] + sellsize[0]/2, cell[1] + sellsize[1]/2
	Up_l = cell[0] - sellsize[0]/2, cell[1] + sellsize[1]/2
	return lowl, lowr, Up_r, Up_l

def points_in_box (pt, bx):
	print (pt, bx)
	for p in pt[:5]:
		# I was working on this function
		print ( p[0])

def main():
	# directory to tiff file
    path_tiff     = "C:/Users/khale/OneDrive/Desktop/Pdal_with_ip/LAS_amsterdam_gf_Laplace.tif"
    ras, cel_Size = read_file(path_tiff)
    # directory to PC file
    path_AHN_PC = "C:/Users/khale/OneDrive/Desktop/Pdal_with_ip/Data_samples/LAS_amsterdam.las"
    AHN_pts = read_PC_Data(path_AHN_PC)
    points_in_box (AHN_pts, ras)
    print (cel_b_box((1,1), (1, 1)))
 

if __name__ == '__main__':
    main()

