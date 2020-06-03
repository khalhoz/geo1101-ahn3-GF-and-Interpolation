### MAIN FILE FOR GROUND-FILTERING/PRE-PROCESSING TESTING ###

# For the documentation, please visit the repo:
# https://github.com/khalhoz/geo1101-ahn3-GF-and-Interpolation

"""You need to specify the folder in which your LAS files are as a command line argument (target folder).
The folder needs to contain the config JSON in the right format (see example file), with name "config.json".
A description of what should be in config.json is found after this block of comments.
The folder also needs to contain a file called "fnames.txt" with the names of the files you want to process,
as also shown in the example file.
Use an ABSOLUTE file path to avoid problems (especially when working from virtual environments)!
Example command line call (in Anaconda Prompt on Windows):
"python C:/Users/geo-geek/some_folder/gf_main.py C:/Users/geo-geek/target_folder/"
The output files will be written to the target folder too. They will be tagged with "_gf".
The returned pdal logs and metadata will also be written there tagged "_log" and "_meta" respectively.

Optionally, you may provide a second argument to specify what to tag the output file with.
For example, if you are running a PDAL pre-rocessing job rather than ground filtering, you could
run "python [target_folder] pre" which would tag the output files with "_pre".
Additionally, a third command line argument may be passed to specify an alternative config file filename.
For example, you can create a separate config for the pre-processin job by calling your file
config_pre.json and passing "python [target_folder] pre pre" in the command line to indicate
that your config file also has the matching "_pre" tag in its filename.
"""

######################################
# DEFAULT JSON PIPELINE CONFIG GUIDE #
######################################
# More info available at https://pdal.io/pipeline.html and also
# in more detail at https://pdal.io/apps/pipeline.html.
# !The default config is based on https://pdal.io/tutorial/ground-filters.html!
# [step 0]  please note we do not have to add a reprojection step
#           to the pipeline - ahn3 is in metres by default, which is
#           what pdal filters need to function properly
# [step 1]  we fill all classifications with 0 (effectively
#           resetting any pre-existing classifications)
# [step 2]  running pdal's noise-elimination algorithm
# [step 3]  running pdal's outlier eleimination algorithm
# [step 4]  running the ground filtering algorithm
#           this uses the algorithm "SMRF" published in Pingel, 2013
#           (we can use filters.pmf for the implementation of Zhang, 2003)
#           parameters: last =  true/false: use only last return or not
#                       ignore = [range]:   which points to let past the
#                                           classification without modification
#                       [...]
#                       (still need to deduce what slope, window,
#                       threshold and scalar are for)
#           EDIT: it appears that "last" has been deprecated? Deleted it from config.
# [step 5]  extracting the point that were classified as ground

from sys import argv
from gf_processing import start_pool

def main():
    if 2 <= len(argv) <= 4: start_pool(True, *argv[1:])
    else: print("Error: Incorrect number of arguments passed. Returning.")

if __name__ == '__main__':
    main()