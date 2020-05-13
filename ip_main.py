### MAIN FILE FOR INTERPOLATION TESTING ###

"""Also uses command line arguments like the PDAL testing script.
The arguments should be specified in this order:
    - target folder (most likely the same as the one you used with PDAL)
    - pixel size (in metres) for interpolation
    - interpolation method, one of:
        - startin-TINlinear
        - startin-Laplace
        - CGAL-NN
        - PDAL-IDW
    - output format, one of:
        - ASC
        - GeoTIFF
    - IDW radius (optional)
    - IDW power (optional)
    - IDW fallback kernel width (optional)
Output files will be dumped in the target folder tagged with the
name of the interpolation method that was used.

Quadrant-based and ellipsodial IDW is not yet implemented,
but radial IDW works well. More info about this in the GitHub readme.
The PDAL-IDW algorithm cannot export in ASC.
"""

from sys import argv
from ip_processing import start_pool

def main():
    if len(argv) == 6 or len(argv) == 5: start_pool(*argv[1:])
    else: print("Error: Incorrect number of arguments passed. Returning.")
    
if __name__ == '__main__':
    main()