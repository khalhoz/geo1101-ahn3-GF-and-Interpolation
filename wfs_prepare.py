### WFS POLYGON IMPORTING AND TILING CODE ###

# For the documentation, please visit the repo:
# https://github.com/khalhoz/geo1101-ahn3-GF-and-Interpolation

import json
from shapely.geometry import shape, Point, LineString, box, Polygon
from shapely.ops import linemerge, unary_union, polygonize
from owslib.wfs import WebFeatureService

def wfs_prepare(bbox, URL):
    """Takes a bounding box and a WFS service URL.
    Requests features in the bounding box, finds polygons that are within
    the bounding box or intersect it. Crops the intersecting geometries
    to the extents of the bounding box, and returns the contained and
    cropped geometries.
    """
    A = Point(bbox[0][0], bbox[1][1])
    B = Point(bbox[0][1], bbox[1][1])
    C = Point(bbox[0][1], bbox[1][0])
    D = Point(bbox[0][0], bbox[1][0])
    bbox_lines = LineString([A,B,C,D,A])
    bbox_object = box(bbox[0][0], bbox[1][0], bbox[0][1], bbox[1][1])
    bag_godzilla = WebFeatureService(url=URL, version='2.0.0')
    response = bag_godzilla.getfeature(typename='BAG3D:pand3d',
                                       bbox=(bbox[0][0], bbox[1][0],
                                             bbox[0][1], bbox[1][1]),
                                       outputFormat='json')
    response_json = json.loads(response.read())
    for feature in response_json['features']:
        rings = feature['geometry']['coordinates']
        for ring_coords, i in zip(rings, range(len(rings))):
            ring = [vx[0:2] for vx in ring_coords]
            feature['geometry']['coordinates'][i] = ring
    out = []
    for feature in response_json['features']:
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
                    feature['geometry']['coordinates'] = [p.exterior.coords]
                    feature['properties']['Shape_Leng'] = p.length
                    feature['properties']['Shape_Area'] = p.area
                    out.append(shape(feature['geometry']))
    return out