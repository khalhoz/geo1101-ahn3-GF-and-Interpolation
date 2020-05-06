### MAIN FILE FOR PDAL TESTING ###

from gf_processing import initialise, start_pool

"""You need to specify the folder in which you LAS files are.
The folder needs to contain the config JSON in the right format (see example file), with name "config.json".
A description of what should be in config.json is found after this block of comments.
The folder also needs to contain a file called "fnames.txt" with the names of the files you want to process,
as also shown in the example file.
Use an ABSOLUTE file path to avoid problems (especially when working from virtual environments)!"""
#####
target_folder = "C:/Users/geo-geek/example_folder/"
#####
"""The output files will be written to the target folder too. They will be tagged with "_out".
# The returned pdal logs and metadata will also be written there tagged "_log" and "_meta" respectively."""

##############################
# JSON PIPELINE CONFIG GUIDE #
##############################
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

def main():
    config, fnames = initialise(target_folder)
    start_pool(target_folder, config, fnames)
    print("Success.")

if __name__ == '__main__':
    main()