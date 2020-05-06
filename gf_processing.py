import json, os, pdal
from time import time
from multiprocessing import Pool, cpu_count

# READ CONFIG
def initialise(target_folder):
    with open(target_folder + "config.json", 'r') as file_in:
        config = file_in.read()
    with open(target_folder + "fnames.txt", 'r') as file_in:
        fnames = file_in.readlines()
    return config, fnames

def worker(mapped):
    """Multiprocessing worker function to be used by the
    p.map function to map objects to, and then start
    multiple times in parallel on separate CPU cores.
    This worker creates and executes the pdal pipeline.
    """
    print("PID {} starting work on {}".format(
        os.getpid(), mapped[1]))
    config, fname = mapped[0], mapped[1]
    config = ('[\n\t"' + fname + '",\n' + config +
              '\n\t"' + fname[:-4] + '_out.las"\n]')
    pipeline = pdal.Pipeline(config)
    start = time()
    pipeline.execute()
    end = time()
    print("PID {} finished processing file.".format(os.getpid()),
          "Time elapsed: {} sec.".format(round(end - start, 2)))
    # TODO:
    # We could write the metadata/logs to disk if we want
    # to keep a better track of what was done to the data.
    #metadata = pipeline.metadata
    #log = pipeline.log  
    # NOTE:
    # The data itself I think could be returned on using:
    #arrays = pipeline.arrays

##############################
# JSON PIPELINE CONFIG GUIDE #
##############################
# !Config based on https://pdal.io/tutorial/ground-filters.html!
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

def start_pool(target_folder, config, fnames):
    """Assembles and executes the worker pool and
    merges the return values (for this test
    code merging is not actually required).
    """
    for i in range(len(fnames)):
        fnames[i] = fnames[i].strip("\n")
        fnames[i] = target_folder + fnames[i]
    processno = cpu_count()
    print("Starting multiprocessing pool on the {}".format(
          processno) + " logical cores found in this PC.")
    if processno < len(fnames):
        print("Warning: more files than processes.\n" +
              "Queue-based multiprocessing not yet implemented, " +
              "ignoring extra files for now.\n")
    elif len(fnames) == 0:
        print("No file names were read. Returning.")
        return
    else: processno = len(fnames)
    pre_map = []
    for i in range(processno): pre_map.append([config, fnames[i]])
    p = Pool(processes = processno)
    out = p.map(worker, pre_map)
    p.close(); p.join()
    print("\nAll workers have returned.")