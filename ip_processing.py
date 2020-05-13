### MULTIPROCESSING POOL-BASED INTERPOLATION CODE ###

import os
from time import time
from multiprocessing import Pool, cpu_count
import numpy as np
from las_prepare import prepare

def execute_startin(pts, res, origin, size, method):
    """Takes the grid parameters and the ground points. Interpolates
    either using the TIN-linear or the Laplace method. Uses a
    -9999 no-data value. Fully based on the startin package.
    """
    import startin
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

def execute_cgal(pts, res, origin, size):
    """Performs CGAL-NN on the input points. Uses something
    that is called regular neighbour coordinates in CGAL, which
    really just means a 2D triangulation with weights.
    The last argument of the query point (qp) is also a weight.
    I don't know what it's for, but it affects the results.
    I set it to zero for now, but we should find out what it is
    during the testing phase.
    """
    from CGAL.CGAL_Kernel import Point_2, Weighted_point_2
    from CGAL.CGAL_Triangulation_2 import Regular_triangulation_2
    from CGAL.CGAL_Interpolation import regular_neighbor_coordinates_2
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

def execute_pdal(target_folder, fpath, size, fmt, rad, pwr, wnd):
    """Sets up a PDAL pipeline that reads a ground filtered LAS
    file, and writes it via GDAL. The GDAL writer has interpolation
    options, exposing the radius, power and a fallback kernel width
    to be configured. More about these in the readme on GitHub.
    """
    import pdal
    if fmt == "ASC":
        print("ASC format for PDAL-IDW is not supported.")
    if fmt == "GeoTIFF":
        config = ('[\n\t"' + fpath + '",\n' +
                  '\n\t{\n\t\t"output_type": "idw"' +
                  ',\n\t\t"resolution": ' + str(size) +
                  ',\n\t\t"radius": ' + str(rad) +
                  ',\n\t\t"power": ' + str(pwr) +
                  ',\n\t\t"window_size": ' + str(wnd) +
                  ',\n\t\t"filename": "' + fpath[:-4] +
                  '_IDW.tif"\n\t}\n]')      
    pipeline = pdal.Pipeline(config)
    pipeline.execute()
    
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

def write_geotiff(raster, origin, fpath):
    """Writes the interpolated TIN-linear and Laplace rasters
    to disk using the GeoTIFF format. The header is based on
    the raster array and a manual definition of the coordinate
    system and an identity affine transform.
    """
    import rasterio
    from rasterio.transform import Affine
    transform = Affine.translation(origin[0], origin[1])
    with rasterio.Env():
        with rasterio.open(fpath, 'w', driver = 'GTiff',
                           height = raster.shape[0],
                           width = raster.shape[1],
                           count = 1,
                           dtype = raster.dtype,
                           crs='EPSG:28992',
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
    size, fpath = mapped[0], (mapped[1] + mapped[2])[:-4] + '_gf.las' 
    start = time()
    if mapped[3] == 'PDAL-IDW':
        execute_pdal(mapped[1], fpath, size, mapped[4],
                     mapped[5], mapped[6], mapped[7])
        end = time()
        print("PID {} finished interpolation and export.".format(os.getpid()),
          "Time elapsed: {} sec.".format(round(end - start, 2)))
        return
    gnd_coords, res, origin = prepare(size, fpath)
    if mapped[3] == 'startin-TINlinear' or mapped[4] == 'startin-Laplace':
        ras = execute_startin(gnd_coords, res, origin, size, mapped[4])
    elif mapped[3] == 'CGAL-NN':
        ras = execute_cgal(gnd_coords, res, origin, size)
    end = time()
    print("PID {} finished interpolation.".format(os.getpid()),
          "Time spent interpolating: {} sec.".format(round(end - start, 2)))
    start = time()
    if mapped[3] == 'startin-TINlinear' and mapped[4] == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_TINlinear.asc')
    if mapped[3] == 'startin-Laplace' and mapped[4] == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_Laplace.asc')
    if mapped[3] == 'CGAL-NN' and mapped[4] == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_NN.asc')
    if mapped[3] == 'startin-TINlinear' and mapped[4] == 'GeoTIFF':
        write_geotiff(ras, origin, fpath[:-4] + '_TINlinear.tif', mapped[3])
    if mapped[3] == 'startin-Laplace' and mapped[4] == 'GeoTIFF':
        write_geotiff(ras, origin, fpath[:-4] + '_Laplace.tif', mapped[3])
    if mapped[3] == 'CGAL-NN' and mapped[4] == 'GeoTIFF':
        write_geotiff(ras, origin, fpath[:-4] + '_NN.tif', mapped[3])
    end = time()
    print("PID {} finished exporting.".format(os.getpid()),
          "Time spent exporting: {} sec.".format(round(end - start, 2)))
    
def start_pool(target_folder, size, method, fmt,
               rad = 10, pwr = 2, wnd = 0):
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
                        method, fmt, rad, pwr, wnd])
    p = Pool(processes = processno)
    p.map(ip_worker, pre_map)
    p.close(); p.join()
    print("\nAll workers have returned.")
    print("Success.")