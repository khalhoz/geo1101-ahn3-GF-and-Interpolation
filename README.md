# Ground Filtering, Interpolation, Hole Filling and Hydro-flattening Testing Environment

## In repo so far:
* `gf_main.py`
* `gf_processing.py`
* `ip_main.py`
* `ip_processing.py`
* target_example (folder)
	* `fnames.txt`
	* `config.json` _(which you can push updates to)_
	* `config_default.json` _(which should not be modified, it is there for reference)_

The testing environment so far includes multiprocessing pool-based implementations of ground filtering via PDAL, TIN-linear and Laplace interpolation via startin, and 
natural neighbour (NN) interpolation via CGAL. I am currently working on integrating a GDAL-based IDW workflow, but it is proving to be more difficult than anticipated.

**Note:** If you are using an Anaconda virtual environment for PDAL/CGAL, you should first active the environment in Anaconda prompt and _then_ run the relevant script
from the same prompt. So, for example:
1. `conda activate [environment_name]`
2. `python [file_path_to_main] [argument_1] [argument_2] [...]`

Relative file paths won't work in virtual environments, so make sure you specify the target folder using a full (absolute) file path.

## PDAL-based ground filtering user guide
1. Copy the example configuration files (`config.json`, `fnames.txt`) to your target folder (in which the LAS files are located).
2. Edit the configuration files.
	1. Specify the names of the LAS files in `fnames.txt` as shown in the example.
	2. Modify the pdal pipeline configuration in `config.json`. Fine-tune the parametrisation to work well with the given data set.
4. Run `gf_main.py` **from the console**. If you run it from an IDE, it will probably not fork the processes properly. It takes one command line argument,
which should be the target folder. An example call in the Windows Anaconda Prompt would be: `python C:/Users/geo-geek/some_folder/gf_main.py C:/Users/geo-geek/target_folder/`

## Interpolation user guide
This is intended to be used after you had generated the ground-filtered files using PDAL (above). The same `fnames.txt` file is used by the program, but it reads the output
files (i.e. files marked `_out.las` rather than the original LAS files).
So, the full workflow is:
1. Check that you still have the PDAL output files and the `fnames.txt` file in your target folder.
2. Run `ip_main.py` **from the console**. If you run it from an IDE, it will probably not fork the processes properly. The following arguments should be provided:
    1. target folder (most likely the same as the one you used with PDAL)
    2. pixel size (in metres) for interpolation
    3. interpolation method, one of:
        * startin-TINlinear
        * startin-Laplace
        * CGAL-NN
    4. output format, one of:
        * ASC
        * GeoTIFF
    5. EPSG override _(optional)_
        * If you wish, you can manually specify an EPSG code for GeoTIFF outputs here, the default is 28992 (Amersfoort).

An example call in the Windows Anaconda Prompt would be: `python C:/Users/geo-geek/some_folder/ip_main.py C:/Users/geo-geek/target_folder/ startin-TINlinear ASC`
The workflow for GDAL-IDW will be slightly different, more info on that later.

### Future work
The current version of this implementation runs as many processes in parallel, as there are input files, hence using all processor cores. This multiprocessing implementation is based
on Python built-in multiprocessing pools. A queue-based implementation would probably work better, but this is something for the scaling group to look at.

Furthermore, we should probably also experiment around with the `filters.pmf` method that is provided by GDAL. The example parametrisation uses `filters.smrf`. You can start experimenting with this simply by changing the JSON parametrisation, the code does not need to be edited.
There's further guidance in `gf_main.py` on what's what in `config.json`.

The GDAL-IDW workflow is currently being worked on, as indicated above. Also relevant to the interpolation is that CGAL-NN takes weights for the query points. I am not sure
conceptually why this is needed and what it is for. It currently uses a weight of zero for each query point, but we should look into what this is for while testing.