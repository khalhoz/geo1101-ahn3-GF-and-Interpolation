### PRE-PROCESSING + INTERPOLATION + POST-PROCESSING SCRIPTS ###

# For the documentation, please visit the repo:
# https://github.com/khalhoz/geo1101-ahn3-GF-and-Interpolation

import os
from time import time
from multiprocessing import Pool, cpu_count
import numpy as np
from las_prepare import las_prepare
from vector_prepare import vector_prepare
from wfs_prepare import wfs_prepare

def triangle_area(a, b, c):
    """Computes the area of a triangle with sides a, b, and c.
    Expects an exception to be raised for problematic area
    calculations, in which case it returns False to indicate
    failure.
    """
    try:
        side1 = np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)
        side2 = np.sqrt((a[0] - c[0])**2 + (a[1] - c[1])**2)
        side3 = np.sqrt((c[0] - b[0])**2 + (c[1] - b[1])**2)
        sp_pa = (side1 + side2 + side3) * 0.5
        return np.sqrt(sp_pa * (sp_pa - side1) *
                       (sp_pa - side2) * (sp_pa - side3))
    except: return False

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
    return ras, tin

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
    s_idx = np.lexsort(pts.T); s_data = pts[s_idx,:]
    mask = np.append([True], np.any(np.diff(s_data[:,:2], axis = 0), 1))
    deduped = s_data[mask]
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

def execute_cgal_cdt(pts, res, origin, size, target_folder):
    """Performs CGAL-CDT on the input points.
    First it removes any potential duplicates from the input points,
    as these would cause issues with the dictionary-based attribute mapping.
    Then, it creates CGAL Point_2 object from these points,
    inserts them into a CGAL Constrained_Delaunay_triangulation_2,
    and then inserts the constraints from shapefiles.
    The constraints do not have elevation values associated.
    Interpolation happens if there are at least two non-constraint
    vertices in a given facet. Otherwise, it either yields the elevation
    of the only non-constraint vertex, or the no-data value is all
    vertices in the facet are constraints.
    It then interpolates (manually, using our code) using TIN-
    linear interpolation via the dictionary-based attribute mapping.
    Extremely long or invalid polygons may mess up the area calculation
    and trigger an exception. These are caught and result in no-data pixels
    which are then filled with values using a median kernel.
    """
    from shapely.geometry import Polygon
    from CGAL.CGAL_Kernel import Point_2
    from CGAL.CGAL_Mesh_2 import Mesh_2_Constrained_Delaunay_triangulation_2
    cdt = Mesh_2_Constrained_Delaunay_triangulation_2()
    s_idx = np.lexsort(pts.T); s_data = pts[s_idx,:]
    mask = np.append([True], np.any(np.diff(s_data[:,:2], axis = 0), 1))
    deduped = s_data[mask]
    cpts = list(map(lambda x: Point_2(*x), deduped[:,:2].tolist()))
    zs = dict(zip([tuple(x) for x in deduped[:,:2]], deduped[:,2]))
    for pt in cpts: cdt.insert(pt)
    poly_fpaths = [
                     'rest_bodies/bbg_rest_of_the_water.shp',
                     'river_bodies/bbg_only_river_bodies.shp',
                     'sea_bodies/bbg_sea_and_big_bodies.shp',
                     # You can add more resources here.
                  ]
    wfs_urls =    [
                     ('http://3dbag.bk.tudelft.nl/data/wfs', 'BAG3D:pand3d'),
                     # You can add more resources here.
                  ]
    in_vecs = []
    for fpath in poly_fpaths:
        vec = vector_prepare([[origin[0], origin[0] + res[0] * size],
                              [origin[1], origin[1] + res[1] * size]],
                             target_folder + fpath)
        if len(vec) != 0: in_vecs.append(vec)
    for wfs in wfs_urls:
        vec = wfs_prepare([[origin[0], origin[0] + res[0] * size],
                           [origin[1], origin[1] + res[1] * size]],
                          wfs[0], wfs[1])
        if len(vec) != 0: in_vecs.append(vec)
    def interpolate(pt):
        tr = cdt.locate(Point_2(pt[0], pt[1]))
        v1 = tr.vertex(0).point().x(), tr.vertex(0).point().y()
        v2 = tr.vertex(1).point().x(), tr.vertex(1).point().y()
        v3 = tr.vertex(2).point().x(), tr.vertex(2).point().y()
        vxs = [v1, v2, v3]
        if (pt[0], pt[1]) in vxs:
            try: zs[(pt[0], pt[1])]
            except: return False
        tr_area = triangle_area(v1, v2, v3)
        if tr_area == False: return False
        ws = [triangle_area((pt[0], pt[1]), v2, v3) / tr_area,
              triangle_area((pt[0], pt[1]), v1, v3) / tr_area,
              triangle_area((pt[0], pt[1]), v2, v1) / tr_area]
        try: vx_zs = [zs[vxs[i]] for i in range(3)]
        except: return False
        return vx_zs[0] * ws[0] + vx_zs[1] * ws[1] + vx_zs[2] * ws[2]
    np.seterr(all='raise')
    for polys in in_vecs:
        for poly in polys:
            if len(poly.exterior.coords[:-1]) < 3: continue
            ring, vals, constraints = [], [], []
            for vx in poly.exterior.coords[:-1]:
                val = interpolate(vx)
                if val == False: continue
                ring.append(vx); vals.append(val)
            try:
                Polygon(ring)
                for val in vals: zs[(vx[0], vx[1])] = val
                for vx in ring:
                    constraints.append(cdt.insert(Point_2(vx[0], vx[1])))
                for vx0, vx1 in zip(constraints, np.roll(constraints, -1)):
                    cdt.insert_constraint(vx0, vx1)
            except: continue
            for interior in poly.interiors:
                ring, vals, constraints = [], [], []
                for vx in interior.coords:
                    val = interpolate(vx)
                    if val == False: continue
                try:
                    Polygon(ring)
                    for val in vals: zs[(vx[0], vx[1])] = val
                    for vx in ring:
                        constraints.append(cdt.insert(Point_2(vx[0], vx[1])))
                    for vx0, vx1 in zip(constraints, np.roll(constraints, -1)):
                        cdt.insert_constraint(vx0, vx1)
                except: continue
    ras = np.zeros([res[1], res[0]])
    yi = 0
    for y in np.arange(origin[1], origin[1] + res[1] * size, size):
        xi = 0
        for x in np.arange(origin[0], origin[0] + res[0] * size, size):
            val = interpolate((x, y))
            if val == False: ras[yi, xi] = -9999
            else: ras[yi, xi] = val
            xi += 1
        yi += 1
    np.seterr(all='warn')
    return ras

def execute_pdal(preproc, target_folder, fpath, size, fmt, rad, pwr, wnd):
    """Sets up a PDAL pipeline that reads a ground filtered LAS
    file, and writes it via GDAL. The GDAL writer has interpolation
    options, exposing the radius, power and a fallback kernel width
    to be configured. More about these in the readme on GitHub.
    """
    import pdal
    if fmt == "GeoTIFF":
        if preproc == False:
            config = ('[\n\t"' + fpath + '",\n' +
                      '\n\t{\n\t\t"output_type": "idw"' +
                      ',\n\t\t"resolution": ' + str(size) +
                      ',\n\t\t"radius": ' + str(rad) +
                      ',\n\t\t"power": ' + str(pwr) +
                      ',\n\t\t"window_size": ' + str(wnd) +
                      ',\n\t\t"filename": "' + fpath[:-4] +
                      '_IDW.tif"\n\t}\n]')
        else:
            with open(target_folder + 
                      "config_preprocess.json", 'r') as file_in:
                preconfig = file_in.read()
            config = ('[\n\t"' + fpath + '",\n' + preconfig +
                      ',\n\t{\n\t\t"output_type": "idw"' +
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
                    ix = tree.query([x, y], rk, tolerance)[1]
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
                rk += incr_rk; i += 1
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
                           dtype = rasterio.float32,
                           crs='EPSG:28992',
                           transform = transform
                           ) as out_file:
            out_file.write(raster.astype(rasterio.float32), 1)

def patch(raster, res, origin, size, min_n):
    """Patches in missing pixel values by applying a median
    kernel (3x3) to estimate its value. This is meant to serve
    as a means of populating missing pixels, not as a means
    of interpolating large areas. The last parameter should
    be an integer that specifies the minimum number of valid
    neighbour values to fill a pixel (0 <= min_n <= 8).
    """
    mp = [[-1, -1], [-1, 0], [-1, 1], [0, -1],
          [0, 1], [1, -1], [1, 0], [1, 1]]
    for yi in range(res[1]):
        for xi in range(res[0]):
            if raster[yi, xi] == -9999:
                vals = []
                for m in range(8):
                    xw, yw = xi + mp[m][0], yi + mp[m][1]
                    if (xw >= 0 and xw < res[0]) and (yw >= 0 and yw < res[1]):
                        val = raster[yw, xw]
                        if val != -9999: vals += [val]
                if len(vals) > min_n: raster[yi, xi] = np.median(vals)

def basic_flattening(target_folder, raster, res, origin, size, tin = False):
    """Reads some pre-determined vector files, tiles them using
    Lisa's code and "burns" them into the output raster. The flat
    elevation of the polygons is estimated by Laplace-interpolating
    at the locations of the polygon vertices. The underlying TIN
    is constructed from the centre points of the raster pixels.
    Rasterisation takes place via rasterio's interface.
    """
    import startin
    from rasterio.features import rasterize
    from rasterio.transform import Affine
    transform = (Affine.translation(origin[0], origin[1])
                 * Affine.scale(size, size))
    x0, x1 = origin[0] + size / 2, origin[0] + ((res[0] - 0.5) * size)
    y0, y1 = origin[1] + size / 2, origin[1] + ((res[1] - 0.5) * size)
    poly_fpaths = [
                     'rest_bodies/bbg_rest_of_the_water.shp',
                     'sea_bodies/bbg_sea_and_big_bodies.shp',
                     # You can add more resources here.
                  ]
    wfs_urls =    [
                     #('http://3dbag.bk.tudelft.nl/data/wfs', 'BAG3D:pand3d'),
                     # You can add more resources here.
                  ]
    in_vecs = []
    for fpath in poly_fpaths:
        vec = vector_prepare([[x0, x1], [y0, y1]], target_folder + fpath)
        if len(vec) != 0: in_vecs.append(vec)
    for wfs in wfs_urls:
        vec = wfs_prepare([[x0, x1], [y0, y1]], wfs[0], wfs[1])
        if len(vec) != 0: in_vecs.append(vec)
    if len(in_vecs) == 0: return
    if tin is False:
        xs, ys = np.linspace(x0, x1, res[0]), np.linspace(y0, y1, res[1])
        xg, yg = np.meshgrid(xs, ys); xg = xg.flatten(); yg = yg.flatten()
        cs = np.vstack((xg, yg, raster.flatten())).transpose()
        data = cs[cs[:,2] != -9999]
        tin = startin.DT(); tin.insert(data)
    elevations = []
    for polys in in_vecs:
        for poly, i in zip(polys, range(len(polys))):
            els = []
            for vx in poly.exterior.coords:
                try: els += [tin.interpolate_laplace(vx[0], vx[1])]
                except: pass
            for interior in poly.interiors:
                for vx in interior.coords:
                    try: els += [tin.interpolate_laplace(vx[0], vx[1])]
                    except: pass
            elevations.append(np.median(els))
    shapes = []
    for polys in in_vecs:
        shapes += [(p, v) for p, v in zip(polys, elevations)]
    raspolys = rasterize(shapes, raster.shape, -9999, transform = transform)
    for yi in range(res[1]):
        for xi in range(res[0]):
            if raspolys[yi, xi] != -9999: raster[yi, xi] = raspolys[yi, xi]
    return tin
            
def hydro_flattening(target_folder, raster, res, origin, size, tin = False):  
    """Reads the river polygons and their skeletons, tiles them using
    Lisa's code and performs hydro-flattening. First cross-sections are
    cast on the skeletons and their elevations are estimated via
    Laplace interpolation. The rivers are then burned into the output
    raster, but each modified pixel is interpolated based on its
    proximity to the closest two cross-sections.
    Since the string of cross-selection is optimised to decrease
    monotonously downstream, under ideal cricumstances this makes the
    interpolation generate a surface with a consistent slope.
    Please read the documentation for more details about the
    algorithm, and its limitations.
    """
    import shapely.geometry as sg
    from rasterio.features import rasterize
    from rasterio.transform import Affine
    transform = (Affine.translation(origin[0], origin[1])
                 * Affine.scale(size, size))
    x0, x1 = origin[0] + size / 2, origin[0] + ((res[0] - 0.5) * size)
    y0, y1 = origin[1] + size / 2, origin[1] + ((res[1] - 0.5) * size)
    river_fpaths = [
                     'river_bodies/bbg_only_river_bodies.shp',
                     # You can add more resources here.
                   ]
    spine_fpaths = [
                     'rivers_final/skeletons_final.shp',
                     # You can add more resources here.
                   ]
    in_rivers = []
    for fpath in river_fpaths:
        vec = vector_prepare([[x0, x1], [y0, y1]], target_folder + fpath)
        if len(vec) != 0: in_rivers.append(vec)
    in_spines = []
    for fpath in spine_fpaths:
        vec = vector_prepare([[x0, x1], [y0, y1]], target_folder + fpath)
        if len(vec) != 0: in_spines.append(vec)
    if len(in_spines) == 0 or len(in_rivers) == 0:
        return
    all_rivers = []
    for river in in_rivers: all_rivers += river
    all_rivers = sg.MultiPolygon(all_rivers)
    cross_sections, elevations, distances = [], [], [0]
    for lstrings in in_spines:
        for lstring in lstrings:
            for i in range(len(lstring.coords[:-1])):
                s0 = lstring.coords[i]; s1 = lstring.coords[i + 1]
                ab = sg.LineString([s0, s1])
                left = ab.parallel_offset(500, 'left')
                right = ab.parallel_offset(500, 'right')
                ortho_r = sg.LineString([s0, right.boundary[1]])
                ortho_l = sg.LineString([s0, left.boundary[0]]) 
                intersections = [all_rivers.boundary.intersection(ortho_r),
                                 all_rivers.boundary.intersection(ortho_l)]
                if type(intersections[0]) == sg.multipoint.MultiPoint:
                    intersections[0] = intersections[0][0]
                elif type(intersections[0]) != sg.point.Point: continue
                if type(intersections[1]) == sg.multipoint.MultiPoint:
                    intersections[1] = intersections[1][0]
                elif type(intersections[1]) != sg.point.Point: continue
                cross_sections.append(sg.LineString([intersections[0],
                                                     intersections[1]]))
                i_vals = []
                try: i_vals.append(tin.interpolate_laplace(s0[0], s0[1]))
                except: pass
                try: i_vals.append(tin.interpolate_laplace(intersections[0].x,
                                                           intersections[0].y))
                except: pass
                try: i_vals.append(tin.interpolate_laplace(intersections[1].x,
                                                           intersections[1].y))
                except: pass
                elevations.append(np.mean(i_vals))
                distances.append(distances[-1] + ab.length)
    if len(elevations) > 5: elevations[0] = np.mean(elevations[0:3])
    elevations = np.array(elevations); distances = np.array(distances[:-1])
    while True:
        mask = np.full(elevations.shape, False); last_valid = elevations[0]
        for i in range(len(elevations) - 1):
            if elevations[i + 1] > last_valid: mask[i + 1] = True
            else: last_valid = elevations[i + 1]
        if True not in mask: break
        elevations[mask] = np.interp(distances[mask], distances[~mask],
                                     elevations[~mask])
    el_dict = {}
    for c_sect, el in zip(cross_sections, elevations):
        el_dict[c_sect.coords[0]] = el
    rshape =  [(all_rivers, 0)]
    raspolys = rasterize(rshape, raster.shape, -9999, transform = transform)
    for yi in range(res[1]):
        for xi in range(res[0]):
            if raspolys[yi, xi] != -9999:
                centre = sg.Point([origin[0] + (xi + 0.5) * size,
                                   origin[1] + (yi + 0.5) * size])
                min_d, closest, next_d, sec_closest = 50000, 0, 50000, 0
                for cs in cross_sections:
                    dist = cs.distance(centre)
                    if dist < min_d:
                        sec_closest = closest; next_d = min_d
                        closest = cs.coords[0]; min_d = dist
                if sec_closest != 0:
                    ele_closest = el_dict[closest]
                    ele_sec_closest = el_dict[sec_closest]
                    u0, w0 = ele_closest, 1 / min_d ** 2
                    u1, w1 = ele_sec_closest, 1 / next_d ** 2
                    asum = u0 * w0 + u1 * w1; bsum = w0 + w1
                    raster[yi, xi] = asum / bsum

def ip_worker(mp):
    """Multiprocessing worker function to be used by the
    p.map function to map objects to, and then start
    multiple times in parallel on separate CPU cores.
    In this case the worker function instances interpolate
    one file each, writing the resulting rasters to disk.
    Runs slightly different workflows depending on the
    desired interpolation method/export format.
    """
    preprocessed, postprocess = mp[0], mp[1]
    size, fpath = mp[2], (mp[3] + mp[4])[:-4] + '_gf.las'
    target_folder, fname, method, fmt = mp[3], mp[4], mp[5], mp[6]
    idw0_polyfpath, idw1, idw2, idw3 = mp[7], mp[8], mp[9], mp[10] 
    idw4, idw5, idw6 = mp[11], mp[12], mp[13]
    print("PID {} starting to work on {}".format(os.getpid(), fname))
    start = time()
    if method == 'PDAL-IDW':
        execute_pdal(preprocessed, target_folder, fpath, size, fmt,
                     idw0_polyfpath, idw1, idw2)
        end = time()
        print("PID {} finished interpolation and export.".format(os.getpid()),
          "Time elapsed: {} sec.".format(round(end - start, 2)))
        return
    if preprocessed == False:
        gnd_coords, res, origin = las_prepare(size, fpath)
    else:
        from math import ceil
        gnd_coords = np.asarray(preprocessed[0].tolist())[:,:3]
        extents = [[min(gnd_coords[:,0]), max(gnd_coords[:,0])],
               [min(gnd_coords[:,1]), max(gnd_coords[:,1])]]
        res = [ceil((extents[0][1] - extents[0][0]) / size),
               ceil((extents[1][1] - extents[1][0]) / size)]
        origin = [np.mean(extents[0]) - (size / 2) * res[0],
                  np.mean(extents[1]) - (size / 2) * res[1]]
    if method == 'startin-TINlinear' or method == 'startin-Laplace':
        ras, tin = execute_startin(gnd_coords, res, origin, size, method)
    elif method == 'CGAL-NN':
        ras = execute_cgal(gnd_coords, res, origin, size)
    elif method == 'CGAL-CDT':
        ras = execute_cgal_cdt(gnd_coords, res, origin, size, target_folder)
    elif method == 'IDWquad':
        ras = execute_idwquad(gnd_coords, res, origin, size,
                              idw0_polyfpath, idw1, idw2, idw3,
                              idw4, idw5, idw6)
    end = time()
    print("PID {} finished interpolating.".format(os.getpid()),
          "Time spent interpolating: {} sec.".format(round(end - start, 2)))
    if postprocess > 0:
        start = time()
        if postprocess == 2 or postprocess == 3 or postprocess == 4:
            if method == 'startin-TINlinear' or method == 'startin-Laplace':
                basic_flattening(target_folder, ras,
                                 res, origin, size, tin)
            else: tin = basic_flattening(target_folder, ras,
                                         res, origin, size)
        if postprocess == 4:
            hydro_flattening(target_folder, ras, res, origin, size, tin)
        if postprocess == 1 or postprocess == 3 or postprocess == 4:
            patch(ras, res, origin, size, 0)
        end = time()
        print("PID {} finished post-processing.".format(os.getpid()),
              "Time spent post-processing: {} sec.".format(
                  round(end - start, 2)))
    start = time()
    if method == 'startin-TINlinear' and fmt == 'GeoTIFF':
        write_geotiff(ras, origin, size, fpath[:-4] + '_TINlinear.tif')
    if method == 'startin-Laplace' and fmt == 'GeoTIFF':
        write_geotiff(ras, origin, size, fpath[:-4] + '_Laplace.tif')
    if method == 'CGAL-NN' and fmt == 'GeoTIFF':
        write_geotiff(ras, origin, size, fpath[:-4] + '_NN.tif')
    if method == 'CGAL-CDT' and fmt == 'GeoTIFF':
        write_geotiff(ras, origin, size, fpath[:-4] + '_TINlinearCDT.tif')
    if method == 'IDWquad' and fmt == 'GeoTIFF':
        write_geotiff(ras, origin, size, fpath[:-4] + '_IDWquad.tif')
    if method == 'startin-TINlinear' and fmt == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_TINlinear.asc')
    if method == 'startin-Laplace' and fmt == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_Laplace.asc')
    if method == 'CGAL-NN' and fmt == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_NN.asc')
    if method == 'CGAL-CDT' and fmt == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_TINlinearCDT.asc')
    if method == 'IDWquad' and fmt == 'ASC':
        write_asc(res, origin, size, ras, fpath[:-4] + '_IDWquad.asc')
    end = time()
    print("PID {} finished exporting.".format(os.getpid()),
          "Time spent exporting: {} sec.".format(round(end - start, 2)))

def start_pool(target_folder, preprocess = "False", postprocess = 0,
               size = 1, method = 'startin-Laplace', fmt = 'GeoTIFF',
               idw0_polyfpath = 5, idw1 = 2, idw2 = 0, idw3 = 2,
               idw4 = 'radial', idw5 = 0.2, idw6 = 3):
    """Assembles and executes the multiprocessing pool.
    The interpolation variants/export formats are handled
    by the worker function (ip_worker(mapped)).
    """
    if preprocess == "False": preprocess = False
    else: preprocess = True
    preprocessed = False
    if preprocess == True and method == 'PDAL-IDW': preprocessed = True
    if int(postprocess) != 0 and method == 'PDAL-IDW':
        print("PDAL-IDW is not yet compatible with post-processing."); return
    elif preprocess == True:
        print("\nRunning pre-processing pool before interpolating.")
        from gf_processing import start_pool as gf_pool
        preprocessed = gf_pool(False, target_folder, '', 'preprocess')
    with open(target_folder + 'fnames.txt', 'r') as file_in:
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
        print("Error: No file names were input. Returning."); return
    pre_map, processno = [], len(fnames)
    if method != 'CGAL-CDT': idw0_polyfpath = float(idw0_polyfpath)
    if preprocess == True and method != 'PDAL-IDW':
        for i in range(processno):
            pre_map.append([preprocessed[i], int(postprocess), float(size),
                            target_folder, fnames[i].strip('\n'), method, fmt,
                            idw0_polyfpath, float(idw1), float(idw2),
                            float(idw3), idw4, float(idw5), float(idw6)])        
    else:
        for i in range(processno):
            pre_map.append([preprocessed, int(postprocess), float(size),
                            target_folder, fnames[i].strip('\n'), method, fmt,
                            idw0_polyfpath, float(idw1), float(idw2),
                            float(idw3), idw4, float(idw5), float(idw6)])
    p = Pool(processes = processno)
    p.map(ip_worker, pre_map)
    p.close(); p.join()
    print("\nAll workers have returned.")