### MAIN FILE FOR INTERPOLATION TESTING ###

"""Also uses command line arguments like the PDAL testing script.
The arguments should be specified in this order:
    - target folder (most likely the same as the one you used with PDAL)
    - pixel size (in metres) for interpolation
    - interpolation method, one of:
        - startin-TINlinear
        - startin-Laplace
        - CGAL-NN
    - output format, one of:
        - ASC
        - GeoTIFF
    - EPSG override
        - If you wish, you can manually specify an EPSG code for GeoTIFF
          outputs here, the default is 28992 (Amersfoort).
Output files will be dumped in the target folder tagged with the
name of the interpolation method that was used.

IDW support is not yet implemented. GDAL is very unsuitable for handling
input in the form of Python/Numpy data structures, so I will need to
come up with an alternative workflow that writes the classified points
to disk in an OGR format rather than LAS, which GDAL can then read,
interpolate and write to disk as GeoTIFF.
"""

from sys import argv
from ip_processing import start_pool

def main():
    if len(argv) == 6 or len(argv) == 5: start_pool(*argv[1:])
    else: print("Error: Incorrect number of arguments passed. Returning.")
    
if __name__ == '__main__':
    main()