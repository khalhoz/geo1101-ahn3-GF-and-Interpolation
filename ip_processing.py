### MULTIPROCESSING POOL-BASED INTERPOLATION CODE ###

import os, math
from time import time
from multiprocessing import Pool, cpu_count
import numpy as np
from laspy.file import File
import rasterio
from rasterio.transform import Affine
import startin
from CGAL.CGAL_Kernel import Point_2, Weighted_point_2
from CGAL.CGAL_Triangulation_2 import Regular_triangulation_2
from CGAL.CGAL_Interpolation import regular_neighbor_coordinates_2

def prepare(size, fpath):
    """Takes the filepath to the input (ground filtered) LAS file, and the
    desired output raster resolution. Reads LAS and outputs the ground points
    as a numpy array. Also establishes some basic raster parameters:
        - the extents
        - the resolution in coordinates
        - the coordinate location of the relative origin (bottom left)
    """
    in_file = File(fpath, mode = "r")
    in_np = np.vstack((in_file.raw_classification,
                       in_file.x, in_file.y, in_file.z)).transpose()
    in_np = in_np[in_np[:,0] == 2].copy()[:,1:]
    extents = [[min(in_np[:,0]), max(in_np[:,0])],
               [min(in_np[:,1]), max(in_np[:,1])]]
    res = [math.ceil((extents[0][1] - extents[0][0]) / size),
           math.ceil((extents[1][1] - extents[1][0]) / size)]
    origin = [np.mean(extents[0]) - (size / 2) * res[0],
              np.mean(extents[1]) - (size / 2) * res[1]]
    return in_np, res, origin

def startin_execute(pts, res, origin, size, method):
    """Takes the grid parameters and the ground points. Simultaneously
    interpolates using the TIN-linear and the Laplace methods. Uses a
    -9999 no-data value. Fully based on the startin package.
    """
    tin = startin.DT(); tin.insert(pts)
    ras = np.zeros(res)
    if method == 'startin-TINlinear':
        def interpolant(x, y): return tin.interpolate_tin_linear(x, y)
    elif method == 'startin-Laplace':
        def interpolant(x, y): return tin.interpolate_laplace(x, y)
    xi = 0
    for x in np.arange(origin[0], origin[0] + res[0] * size, size):
        yi = 0
        for y in np.arange(origin[1] + res[1] * size, origin[1], -size):
            tri = tin.locate(x, y)
            if tri != [] and 0 not in tri:
                ras[xi, yi] = interpolant(x, y)
            else: ras[xi, yi] = -9999
            yi += 1
        xi += 1
    return ras

def cgal_execute(pts, res, origin, size):
    """Performs CGAL-NN on the input points. Uses something
    that is called regular neighbour coordinates in CGAL, which
    really just means a 2D triangulation with weights.
    The last argument of the query point (qp) is also a weight.
    I don't know what it's for, but it affects the results.
    I set it to zero for now, but we should find out what it is
    during the testing phase.
    """
    cpts = list(map(lambda x: Point_2(*x), pts[:,:2].tolist()))
    pairs = [(p, z) for p, z in zip(cpts, pts[:,2].tolist())]
    wpts = map(lambda x: Weighted_point_2(*x), pairs)
    tin = Regular_triangulation_2()
    for pt in wpts: tin.insert(pt)
    ras = np.zeros(res)
    xi = 0
    for x in np.arange(origin[0], origin[0] + res[0] * size, size):
        yi = 0
        for y in np.arange(origin[1] + res[1] * size, origin[1], -size):
            qp = Weighted_point_2(Point_2(x, y), 0)
            interpvalue = regular_neighbor_coordinates_2(tin, qp, [])
            if interpvalue[1] == True: ras[xi, yi] = interpvalue[0]
            else: ras[xi, yi] = -9999
            yi += 1
        xi += 1
    return ras

def write_asc(res, origin, size, raster, fpath):
    """Writes the interpolated TIN-linear and Laplace rasters
    to disk using the ASC format. The header is based on the
    pre-computed raster parameters.
    """
    with open(fpath, "w") as file_out:
        file_out.write("NCOLS " + str(res[0]) + "\n")
        file_out.write("NROWS " + str(res[1]) + "\n")
        file_out.write("XLLCORNER " + str(origin[0]) + "\n")
        file_out.write("YLLCORNER " + str(origin[1]) + "\n")
        file_out.write("CELLSIZE " + str(size) + "\n")
        file_out.write("NODATA_VALUE " + str(-9999) + "\n")
        for yi in range(res[1]):
            for xi in range(res[0]):
                file_out.write(str(raster[xi, yi]) + " ")
            file_out.write("\n")

def write_geotiff(raster, origin, fpath, epsg):
    """Writes the interpolated TIN-linear and Laplace rasters
    to disk using the GeoTIFF format. The header is based on
    the raster array and a manual definition of the coordinate
    system and an identity affine transform.
    """
    transform = Affine.translation(origin[0], origin[1])
    with rasterio.Env():
        with rasterio.open(fpath, 'w', driver = 'GTiff',
                           height = raster.shape[0],
                           width = raster.shape[1],
                           count = 1,
                           dtype = raster.dtype,
                           crs='EPSG:' + epsg,
                           transform = transform
                           ) as out_file:
            out_file.write(raster, 1)
           
def ip_worker(mapped):
    """Multiprocessing worker function to be used by the
    p.map function to map objects to, and then start
    multiple times in parallel on separate CPU cores.
    In this case the worker function instances interpolate
    one file each, writing the resulting rasters to disk.
    Runs slightly different workflows depending on the
    desired interpolation method/export format.
    """
    print("PID {} starting to interpolate file {}".format(
        os.getpid(), mapped[2]))
    size, fpath = mapped[0], (mapped[1] + mapped[2])[:-4] + '_out.las' 
    start = time()
    gnd_coords, res, origin = prepare(size, fpath)
    if mapped[4] == 'startin-TINlinear' or mapped[4] == 'startin-Laplace':
        ras = startin_execute(gnd_coords, res, origin, size, mapped[4])
    elif mapped[4] == 'CGAL-NN':
        ras = cgal_execute(gnd_coords, res, origin, size)
    elif mapped[4] == 'GDAL-IDW-radial':
        pass
    end = time()
    print("PID {} finished interpolation.".format(os.getpid()),
          "Time spent interpolating: {} sec.".format(round(end - start, 2)))
    start = time()
    if mapped[4] == 'startin-TINlinear' and mapped[5] == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_TINlinear.asc')
    if mapped[4] == 'startin-Laplace' and mapped[5] == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_Laplace.asc')
    if mapped[4] == 'startin-TINlinear' and mapped[5] == 'GeoTIFF':
        write_geotiff(ras, origin, fpath[:-4] + '_TINlinear.tif', mapped[3])
    if mapped[4] == 'startin-Laplace' and mapped[5] == 'GeoTIFF':
        write_geotiff(ras, origin, fpath[:-4] + '_Laplace.tif', mapped[3])
    if mapped[4] == 'CGAL-NN' and mapped[5] == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_NN.asc')
    if mapped[4] == 'CGAL-NN' and mapped[5] == 'GeoTIFF':
        write_geotiff(ras, origin, fpath[:-4] + '_NN.tif', mapped[3])
    end = time()
    print("PID {} finished exporting.".format(os.getpid()),
          "Time spent exporting: {} sec.".format(round(end - start, 2)))
    
def start_pool(target_folder, size, method, fmt, epsg = '28992'):
    """Assembles and executes the multiprocessing pool.
    The interpolation variants/export formats are handled
    by the worker function (ip_worker(mapped)).
    """
    with open(target_folder + "fnames.txt", 'r') as file_in:
        fnames = file_in.readlines()
    cores = cpu_count()
    print("\nStarting interpolation pool of processes on the {}".format(
        cores) + " logical cores found in this PC.\n")
    if cores < len(fnames):
        print("Warning: more processes in pool than processor cores.\n" +
              "Optimally, roughly as many processes as processor " +
              "cores should be run concurrently.\nYou are starting " +
              str(len(fnames)) + " processes on " + str(cores) + " cores.\n")
    elif len(fnames) == 0:
        print("Error: No file names were read. Returning."); return
    processno = len(fnames)
    pre_map = []
    for i in range(processno):
        pre_map.append([float(size), target_folder, fnames[i].strip("\n"),
                        epsg, method, fmt])
    p = Pool(processes = processno)
    p.map(ip_worker, pre_map)
    p.close(); p.join()
    print("\nAll workers have returned.")
    print("\nSuccess.")