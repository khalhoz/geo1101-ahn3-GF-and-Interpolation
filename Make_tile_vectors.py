from shapely.geometry import Point, LineString, box, Polygon
from shapely.ops import linemerge, unary_union, polygonize
import laspy
import fiona
import os

def split_the_tiles(lasfile_path,name):
    #get the bounding box
    inFile = laspy.file.File(lasfile_path, mode="r")
    header = inFile.header
    min_x = header.min[0]
    max_x = header.max[0]
    min_y = header.min[1]
    max_y = header.max[1]

    #make the 4 corner points and their "cut" line
    A = Point(min_x, max_y)
    B = Point(max_x, max_y)
    C = Point(max_x, min_y)
    D = Point(min_x, min_y)
    boundary_line = LineString([A,B,C,D,A])
    bbox = box(min_x,min_y,max_x,max_y)


    def get_tile(box,boundary_line, input, output):
        for feature in input:
            found_intersect = False
            merge_list= [boundary_line]
            coord = feature['geometry']['coordinates']
            for ring in coord:
                poly = Polygon(ring)
                if poly.within(box):
                    output.write(feature)
                elif poly.intersects(box):
                    found_intersect = True
                    merge_list.append(poly.boundary)



            if found_intersect:
                if len(coord) > 1:
                    poly = Polygon(coord[0], coord[1:])
                else:
                    poly = Polygon(coord[0])
                merged = linemerge(merge_list)
                borders = unary_union(merged)
                polygons = polygonize(borders)
                for p in polygons:
                    if p.within(box) and poly.contains(p.buffer(-1e-8)):
                        print("TRUE")
                        feature['geometry']['coordinates']= [p.exterior.coords]
                        feature['properties']['Shape_Leng'] = p.length
                        feature['properties']['Shape_Area'] = p.area
                        output.write(feature)

    if not os.path.exists(name):
        os.makedirs(name)
        parentdir = name
        path_river = os.path.join(parentdir, "river_bodies")
        os.makedirs(path_river)
        path_sea = os.path.join(parentdir, "sea_bodies")
        os.makedirs(path_sea)
        path_rest = os.path.join(parentdir, "rest_bodies")
        os.makedirs(path_rest)

    with fiona.open('river_bodies/bbg_only_river_bodies.shp') as src_river:
        # we safe the metadata as parameter, seeing we want to re-use those in the sub files
        meta = src_river.meta
        # Open the output files, using the same format driver and coordinate reference system as the source (**meta).
        with fiona.open(str(name) +'/river_bodies/tile_river_bodies.shp', 'w', **meta) as river:
            get_tile(bbox,boundary_line,src_river,river)

    with fiona.open('sea_bodies/bbg_sea_and_big_bodies.shp') as src_sea:
        # we safe the metadata as parameter, seeing we want to re-use those in the sub files
        meta = src_sea.meta
        # Open the output files, using the same format driver and coordinate reference system as the source (**meta).
        with fiona.open(str(name)+'/sea_bodies/tile_sea_bodies.shp', 'w', **meta) as sea:
            get_tile(bbox,boundary_line,src_sea,sea)

    with fiona.open('rest_bodies/bbg_rest_of_the_water.shp') as src:
        # we safe the metadata as parameter, seeing we want to re-use those in the sub files
        meta = src.meta
        # Open the output files, using the same format driver and coordinate reference system as the source (**meta).
        with fiona.open(str(name)+'/rest_bodies/tile_rest_bodies.shp', 'w', **meta) as rest:
            get_tile(bbox,boundary_line,src,rest)


split_the_tiles("test_files_las/5 samples/LAS_biesbosh.las", "biesbosch")