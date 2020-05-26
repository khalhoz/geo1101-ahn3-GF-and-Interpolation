### MULTIPROCESSING POOL-BASED PDAL GROUND FILTERING CODE ###
import numpy as np
import os
from time import time
from multiprocessing import Pool, cpu_count
import pdal

def initialise(target_folder, json):
    """This function reads and returns the configuration JSON file and the
    list of file names. These are assumed to be called config.json
    and fnames.txt and to be found in the same folder as the LAS files.
    """
    with open(target_folder + "config" + json + ".json", 'r') as file_in:
        config = file_in.read()
    with open(target_folder + "fnames.txt", 'r') as file_in:
        fnames = file_in.readlines()
    return config, fnames

def worker(mapped):
    """Multiprocessing worker function to be used by the
    p.map function to map objects to, and then start
    multiple times in parallel on separate CPU cores.
    In this case the worker function instances ground
    filter one file each, and return the resulting log,
    metadata and the point cloud data array itself.
    """
    print("PID {} starting to ground filter file {}".format(
        os.getpid(), mapped[2]))
    config, fpath, write = mapped[0], mapped[1] + mapped[2], mapped[3]
    if write == True:
        tag = mapped[4]
        config = ('[\n\t"' + fpath + '",\n' + config +
                  '\n\t"' + fpath[:-4] + '_' + tag + '.las"\n]')
    else: config = ('[\n\t"' + fpath + '",\n' + config + '\n]')
    pipeline = pdal.Pipeline(config)
    start = time()
    pipeline.execute()
    end = time()
    print("PID {} finished ground filtering.".format(os.getpid()),
          "Time elapsed: {} sec.".format(round(end - start, 2)))
    log = pipeline.log
    metadata = pipeline.metadata
    arrays = pipeline.arrays
    return log, metadata, arrays

def start_pool(write, target_folder, tag = "gf", json = ""):
    """Assembles and executes the multiprocessing pool and
    merges the return values (logs, metadata, data arrays).
    It currently writes the logs and metadata to disk and
    does not actually do anything with the returned data arrays.
    As we discussed with Maarten on GitHub, maybe it would be
    better to pass on the ground filtered data to interpolation
    intelligently rather than to write all the ground filtered LAS
    files to disk (as it is done by this script currently).
    In that case the list array_gf would need to be passed on,
    as it collects all the returned output data arrays of all
    the ground filtering processes completed as part of the
    multiprocessing pool.
    """
    if json != "": json = "_" + json
    config, fnames = initialise(target_folder, json)
    cores = cpu_count()
    print("\nStarting GF/preprocessing pool of processes on the {}".format(
        cores) + " logical cores found in this PC.\n")
    if cores < len(fnames):
        print("Warning: more processes in pool than processor cores.\n" +
              "Optimally, roughly as many processes as processor " +
              "cores should be run concurrently.\nYou are starting " +
              str(len(fnames)) + " processes on " + str(cores) + " cores.\n")
    elif len(fnames) == 0:
        print("Error: No file names were read. Returning."); return
    processno = len(fnames)
    pre_map = []
    for i in range(processno):
        fnames[i] = fnames[i].strip("\n")
        pre_map.append([config, target_folder, fnames[i], bool(write), tag])
    p = Pool(processes = processno)
    out = p.map(worker, pre_map)
    logs_gf, meta_gf, arrays_gf = [], [], []
    for returned in out:
        logs_gf.append(returned[0])
        meta_gf.append(returned[1])
        arrays_gf.append(returned[2])
    p.close(); p.join()
    print("\nAll workers have returned. Writing logs and metadata.")
    # The logs appear to be empty. I don't know whether this is
    # because the current version of pdal does not use these logs,
    # or because something needs to be configured manually.
    for fname, log, meta in zip(fnames, logs_gf, meta_gf):
        with open(target_folder + fname[:-4] + "_meta.json", "w") as file_out:
            file_out.write(meta)
        with open(target_folder + fname[:-4] + "_log.txt", "w") as file_out:
            file_out.write(log)
    return arrays_gf
    print("Success.")