### MAIN FILE FOR INTERPOLATION TESTING ###

"""Also uses command line arguments like the PDAL testing script.
The arguments should be specified in this order:
    - target folder (most likely the same as the one you used with PDAL)
      [compulsory, no default value]
    - pixel size (in metres) for interpolation
      [default: 1 metre]
    - interpolation method, one of:
        - startin-TINlinear
        - [default] startin-Laplace
        - CGAL-NN
        - PDAL-IDW
        - IDWquad
    - output format, one of:
        - ASC
        - [default] GeoTIFF
    - IDW radius (for PDAL-IDW)
      // STARTING IDW radius/number of neighbours to query (for IDWquad)
    - IDW power (for PDAL-IDW and IDWquad)
    - IDW fallback kernel width (for PDAL-IDW)
      // MINIMUM number of points per quadrant (for IDWquad)
    - radius/number of neighbours INCREMENT value (for IDWquad)
    - IDWquad method, one of:
        - radial
        - k-nearest
    - IDWquad tolerance (epsilon)
    - IDWquad maximum number of iteration before declaring no-data

All IDW parameters are optional, but it is assumed the user will fine-
tune them, hence the defaults are not listed.

Output files will be dumped in the target folder tagged with the
name of the interpolation method that was used.
"""

from sys import argv
from ip_processing import start_pool

def main():
    if 3 <= len(argv) <= 10: start_pool(*argv[1:])
    else: print("Error: Incorrect number of arguments passed. Returning.")
    
if __name__ == '__main__':
    main()