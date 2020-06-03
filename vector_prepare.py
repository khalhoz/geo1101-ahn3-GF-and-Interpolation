### VECTOR FILE IMPORTING AND VECTOR TILING CODE ###

from shapely.geometry import shape, Point, LineString, box, Polygon
from shapely.ops import linemerge, unary_union, polygonize
import fiona

def vector_prepare(bbox, fpath):
    A = Point(bbox[0][0], bbox[1][1])
    B = Point(bbox[0][1], bbox[1][1])
    C = Point(bbox[0][1], bbox[1][0])
    D = Point(bbox[0][0], bbox[1][0])
    boundary_line = LineString([A,B,C,D,A])
    bbox_object = box(bbox[0][0], bbox[1][0], bbox[0][1], bbox[1][1])
    out = []
    for feature in fiona.open(fpath):
        merger = [boundary_line]
        coords = feature['geometry']['coordinates']
        for ring in coords:
            poly = Polygon(ring)
            if poly.within(bbox_object):
                out.append(shape(feature['geometry']))
            elif poly.intersects(bbox_object):
                merger.append(poly.boundary)
        if len(merger) != 1:
            if len(coords) > 1: poly = Polygon(coords[0], coords[1:])
            else: poly = Polygon(coords[0])
            merged = linemerge(merger)
            borders = unary_union(merged)
            polygons = polygonize(borders)
            for p in polygons:
                if p.within(bbox_object) and poly.contains(p.buffer(-1e-8)):
                    feature['geometry']['coordinates']= [p.exterior.coords]
                    feature['properties']['Shape_Leng'] = p.length
                    feature['properties']['Shape_Area'] = p.area
                    out.append(shape(feature['geometry']))
    return out