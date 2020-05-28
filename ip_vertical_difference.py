import rasterio
from rasterio.transform import Affine
import laspy 
import numpy as np
from multiprocessing import Pool


def read_file (directory):
    src = rasterio.open(directory)
    width = src.width
    height = src.height 
    bound = src.bounds
    raster_array = src.read(1)
    cell_size = ((src.transform * (src.width, src.height)) [0] - (src.transform * (0, 0))[0])/width, ((src.transform * (src.width, src.height)) [1] - (src.transform * (0, 0))[1])/height
    return [raster_array, cell_size, bound, width, height]


def read_PC_Data (directory): 
    pc = laspy.file.File(directory, mode='r')
    pts = np.vstack((pc.x,pc.y,pc.z)).transpose()
    return pts


def calculate_differences(pc_pts, raster_array, raster_cell_size, raster_bbox, raster_width, raster_height):
    verticle_differences = np.zeros((raster_height, raster_width))
    xmin = raster_bbox[0]
    ymin = raster_bbox[1]
    sizex = raster_cell_size[0]
    sizey = raster_cell_size[1]

    for pt in pc_pts:
        col = int((pt[0]-xmin)/sizex)
        row = int((pt[1]-ymin)/sizey)
        if raster_array[row][col] == -9999: z = pt[2]
        else: z = pt[2]-raster_array[row][col]
        verticle_differences[row][col] += z
        
    for col in range(raster_width):
        for row in range(raster_height):
            if verticle_differences[row][col] == 0 and raster_array[row][col] != -9999:
                verticle_differences[row][col] = raster_array[row][col]

    return verticle_differences


def main(AHN_pc_file):
    AHN_pc = read_PC_Data("AHN/" + AHN_pc_file + ".las")
    interpolation_methods = ["IDW", "NN", "Laplace", "TINlinear"]
    file_name = "DSM_" + AHN_pc_file[4:] + "_"
    for method in interpolation_methods:
        raster = method + "/" + file_name + method + ".tif"
        raster_info = read_file(raster)
        differences_values = calculate_differences(AHN_pc, raster_info[0], raster_info[1], raster_info[2], raster_info[3], raster_info[4])
        transform = (Affine.translation(raster_info[2][0], raster_info[2][3]) * Affine.scale(raster_info[1][0], raster_info[1][1]))
        with rasterio.Env():
            with rasterio.open(method + "/" + file_name +'_verticle_differences.tif', 'w', driver = 'GTiff',
                        height = raster_info[4], width = raster_info[3], count = 1,
                        dtype = differences_values.dtype, crs='EPSG:28992', transform = transform) as dst:
                dst.write(differences_values, 1)


if __name__ == '__main__':
    pre_map = ["LAS_amsterdam", "LAS_biesbosh", "LAS_Delft", "LAS_Groningen", "LAS_Veluwezoom"]
    p = Pool(processes = 5)
    p.map(main, pre_map)
    p.close(); p.join() 
    print("Done")