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
    ras = np.zeros([res[1], res[0]])
    if method == 'startin-TINlinear':
        def interpolant(x, y): return tin.interpolate_tin_linear(x, y)
    elif method == 'startin-Laplace':
        def interpolant(x, y): return tin.interpolate_laplace(x, y)
    yi = 0
    for y in np.arange(origin[1], origin[1] + res[1] * size, size):
        xi = 0
        for x in np.arange(origin[0], origin[0] + res[0] * size, size):
            tri = tin.locate(x, y)
            if tri != [] and 0 not in tri:
                ras[yi, xi] = interpolant(x, y)
            else: ras[yi, xi] = -9999
            xi += 1
        yi += 1
    return ras

def execute_cgal(pts, res, origin, size):
    """Performs CGAL-NN on the input points.
    First it removes any potential duplicates from the
    input points, as these would cause issues with the
    dictionary-based attribute mapping.
    Then, it creates CGAL Point_2 object from these points,
    inserts them into a CGAL Delaunay_triangulation_2, and
    performs interpolation using CGAL natural_neighbor_coordinate_2
    by finding the attributes (Z coordinates) via the dictionary
    that was created from the deduplicated points.
    """
    from CGAL.CGAL_Kernel import Point_2
    from CGAL.CGAL_Triangulation_2 import Delaunay_triangulation_2
    from CGAL.CGAL_Interpolation import natural_neighbor_coordinates_2
    sorted_idx = np.lexsort(pts.T)
    sorted_data = pts[sorted_idx,:]
    row_mask = np.append([True], np.any(np.diff(sorted_data[:,:2],
                                                axis = 0), 1))
    deduped = sorted_data[row_mask]
    cpts = list(map(lambda x: Point_2(*x), deduped[:,:2].tolist()))
    zs = dict(zip([tuple(x) for x in deduped[:,:2]], deduped[:,2]))
    tin = Delaunay_triangulation_2()
    for pt in cpts: tin.insert(pt)
    ras = np.zeros([res[1], res[0]])
    yi = 0
    for y in np.arange(origin[1], origin[1] + res[1] * size, size):
        xi = 0
        for x in np.arange(origin[0], origin[0] + res[0] * size, size):
            nbrs = [];
            qry = natural_neighbor_coordinates_2(tin, Point_2(x, y), nbrs)
            if qry[1] == True:
                z_out = 0
                for nbr in nbrs:
                    z, w = zs[(nbr[0].x(), nbr[0].y())], nbr[1] / qry[0]
                    z_out += z * w
                ras[yi, xi] = z_out
            else: ras[yi, xi] = -9999
            xi += 1
        yi += 1
    return ras

def execute_pdal(target_folder, fpath, size, fmt, rad, pwr, wnd):
    """Sets up a PDAL pipeline that reads a ground filtered LAS
    file, and writes it via GDAL. The GDAL writer has interpolation
    options, exposing the radius, power and a fallback kernel width
    to be configured. More about these in the readme on GitHub.
    """
    import pdal
    if fmt == "GeoTIFF":
        config = ('[\n\t"' + fpath + '",\n' +
                  '\n\t{\n\t\t"output_type": "idw"' +
                  ',\n\t\t"resolution": ' + str(size) +
                  ',\n\t\t"radius": ' + str(rad) +
                  ',\n\t\t"power": ' + str(pwr) +
                  ',\n\t\t"window_size": ' + str(wnd) +
                  ',\n\t\t"filename": "' + fpath[:-4] +
                  '_IDW.tif"\n\t}\n]')      
        pipeline = pdal.Pipeline(config); pipeline.execute()
    elif fmt == "ASC": print("ASC format for PDAL-IDW is not supported.")

def execute_idwquad(pts, res, origin, size,
                    start_rk, pwr, minp, incr_rk, method, tolerance, maxiter):
    """Creates a KD-tree representation of the tile's points and
    executes a quadrant-based IDW algorithm on them. Although the
    KD-tree is based on a C implementation, the rest is coded in
    pure Python (below). Keep in mind that because of this, this
    is inevitably slower than the rest of the algorithms here.
    To optimise performance, one is advised to fine-tune the
    parametrisation, especially tolerance and maxiter.
    More info in the GitHub readme.
    """
    from scipy.spatial import cKDTree
    ras = np.zeros([res[1], res[0]])
    tree = cKDTree(np.array([pts[:,0], pts[:,1]]).transpose())
    yi = 0
    for y in np.arange(origin[1], origin[1] + res[1] * size, size):
        xi = 0
        for x in np.arange(origin[0], origin[0] + res[0] * size, size):
            done, i, rk = False, 0, start_rk
            while done == False:
                if method == "radial":
                    ix = tree.query_ball_point([x, y], rk, tolerance)
                elif method == "k-nearest":
                    ix = tree.query([x, y], rk, tolerance)
                xyp = pts[ix]
                qs = [
                        xyp[(xyp[:,0] < x) & (xyp[:,1] < y)],
                        xyp[(xyp[:,0] > x) & (xyp[:,1] < y)],
                        xyp[(xyp[:,0] < x) & (xyp[:,1] > y)],
                        xyp[(xyp[:,0] > x) & (xyp[:,1] > y)]
                     ]
                if min(qs[0].size, qs[1].size,
                       qs[2].size, qs[3].size) >= minp: done = True
                elif i == maxiter:
                    ras[yi, xi] = -9999; break
                rk += incr_rk
                i += 1
            else:
                asum, bsum = 0, 0
                for pt in xyp:
                    dst = np.sqrt((x - pt[0])**2 + (y - pt[1])**2)
                    u, w = pt[2], 1 / dst ** pwr
                    asum += u * w; bsum += w
                    ras[yi, xi] = asum / bsum
            xi += 1
        yi += 1
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
        for yi in range(res[1] - 1, -1, -1):
            for xi in range(res[0]):
                file_out.write(str(raster[yi, xi]) + " ")
            file_out.write("\n")

def write_geotiff(raster, origin, size, fpath):
    """Writes the interpolated TIN-linear and Laplace rasters
    to disk using the GeoTIFF format. The header is based on
    the raster array and a manual definition of the coordinate
    system and an identity affine transform.
    """
    import rasterio
    from rasterio.transform import Affine
    transform = (Affine.translation(origin[0], origin[1])
                 * Affine.scale(size, size))
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

def ip_worker(mp):
    """Multiprocessing worker function to be used by the
    p.map function to map objects to, and then start
    multiple times in parallel on separate CPU cores.
    In this case the worker function instances interpolate
    one file each, writing the resulting rasters to disk.
    Runs slightly different workflows depending on the
    desired interpolation method/export format.
    """
    size, fpath, fname = mp[0], (mp[1] + mp[2])[:-4] + '_gf.las', mp[2]
    target_folder, method, fmt = mp[1], mp[3], mp[4]
    idw0, idw1, idw2, idw3 = mp[5], mp[6], mp[7], mp[8] 
    idw4, idw5, idw6 = mp[9], mp[10], mp[11]
    print("PID {} starting to interpolate file {}".format(os.getpid(), fname))
    start = time()
    if method == 'PDAL-IDW':
        execute_pdal(target_folder, fpath, size, fmt, idw0, idw1, idw2)
        end = time()
        print("PID {} finished interpolation and export.".format(os.getpid()),
          "Time elapsed: {} sec.".format(round(end - start, 2)))
        return
    gnd_coords, res, origin = prepare(size, fpath)
    if method == 'startin-TINlinear' or method == 'startin-Laplace':
        ras = execute_startin(gnd_coords, res, origin, size, method)
    elif method == 'CGAL-NN':
        ras = execute_cgal(gnd_coords, res, origin, size)
    elif method == 'IDWquad':
        ras = execute_idwquad(gnd_coords, res, origin, size,
                              idw0, idw1, idw2, idw3, idw4, idw5, idw6)
    end = time()
    print("PID {} finished interpolation.".format(os.getpid()),
          "Time spent interpolating: {} sec.".format(round(end - start, 2)))
    start = time()
    if method == 'startin-TINlinear' and fmt == 'GeoTIFF':
        write_geotiff(ras, origin, size, fpath[:-4] + '_TINlinear.tif')
    if method == 'startin-Laplace' and fmt == 'GeoTIFF':
        write_geotiff(ras, origin, size, fpath[:-4] + '_Laplace.tif')
    if method == 'CGAL-NN' and fmt == 'GeoTIFF':
        write_geotiff(ras, origin, size, fpath[:-4] + '_NN.tif')
    if method == 'IDWquad' and fmt == 'GeoTIFF':
        write_geotiff(ras, origin, size, fpath[:-4] + '_IDWquad.tif')
    if method == 'startin-TINlinear' and fmt == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_TINlinear.asc')
    if method == 'startin-Laplace' and fmt == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_Laplace.asc')
    if method == 'CGAL-NN' and fmt == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_NN.asc')
    if method == 'IDWquad' and fmt == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_IDWquad.asc')
    end = time()
    print("PID {} finished exporting.".format(os.getpid()),
          "Time spent exporting: {} sec.".format(round(end - start, 2)))

def start_pool(target_folder, size = 1, method = "startin-Laplace",
               fmt = "GeoTIFF", idw0 = 5, idw1 = 2, idw2 = 0,
               idw3 = 2, idw4 = "radial", idw5 = 0.2, idw6 = 3):
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
                        method, fmt, float(idw0), float(idw1), float(idw2),
                        float(idw3), idw4, float(idw5), float(idw6)])
    p = Pool(processes = processno)
    p.map(ip_worker, pre_map)
    p.close(); p.join()
    print("\nAll workers have returned.")
    print("Success.")