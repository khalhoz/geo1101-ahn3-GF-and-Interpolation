### VECTOR FILE IMPORTING AND VECTOR TILING CODE ###

# For the documentation, please visit the repo:
# https://github.com/khalhoz/geo1101-ahn3-GF-and-Interpolation

from shapely.geometry import Point, LineString, Polygon, shape, box
from shapely.ops import split, linemerge, unary_union, polygonize
import fiona

def vector_prepare(bbox, fpath):
    """Takes a bounding box and a file path to a vector file.
    Reads the vector file, finds polygons that are within the
    bounding box or intersect it. Crops the intersecting geometries
    to the extents of the bounding box, and returns the contained and
    cropped geometries.
    """
    A = Point(bbox[0][0], bbox[1][1])
    B = Point(bbox[0][1], bbox[1][1])
    C = Point(bbox[0][1], bbox[1][0])
    D = Point(bbox[0][0], bbox[1][0])
    bbox_lines = LineString([A,B,C,D,A])
    bbox_object = box(bbox[0][0], bbox[1][0], bbox[0][1], bbox[1][1])
    out = []
    for feature in fiona.open(fpath):
        if feature['geometry']['type'] == 'Polygon':
            merger = [bbox_lines]
            rings = feature['geometry']['coordinates']
            for ring_coords in rings:
                ring = Polygon(ring_coords)
                if ring.within(bbox_object): out.append(shape(feature['geometry']))
                elif ring.intersects(bbox_object): merger.append(ring.boundary)
            if len(merger) != 1:
                if len(rings) > 1: poly = Polygon(rings[0], rings[1:])
                else: poly = Polygon(rings[0])
                merged = linemerge(merger)
                borders = unary_union(merged)
                polygons = polygonize(borders)
                for p in polygons:
                    if p.within(bbox_object) and poly.contains(p.buffer(-1e-8)):
                        feature['geometry']['coordinates']= [p.exterior.coords]
                        feature['properties']['Shape_Leng'] = p.length
                        feature['properties']['Shape_Area'] = p.area
                        out.append(shape(feature['geometry']))
        elif feature['geometry']['type'] == 'LineString':
            lstr = LineString(feature['geometry']['coordinates'])
            if lstr.within(bbox_object): out.append(shape(lstr))
            elif lstr.intersects(bbox_object):
                lstrs = split(lstr, bbox_object)
                for l in lstrs:
                    if l.within(bbox_object):
                        feature['geometry']['coordinates']= l.coords
                        feature['properties']['Shape_Leng'] = l.length
                        feature['properties']['Shape_Area'] = l.area
                        out.append(shape(feature['geometry']))
    return out