# Ground Filtering, Interpolation, Hole Filling and Hydro-flattening Testing Environment

**This is the Python repo of the ground filtering/interpolation team of the AHN3 GEO1101 (synthesis project) group of 2020.**

## In repo so far:
* `las_prepare.py` _(factored out from `gf_processing.py`)_
* `gf_main.py`
* `gf_processing.py`
* `ip_main.py`
* `ip_processing.py`
* `gdal_attempt.py`
* target_example (folder)
	* `fnames.txt`
	* `config.json` _(which you can push updates to)_
	* `config_default.json` _(which should not be modified, it is there for reference)_

The testing environment so far includes multiprocessing pool-based implementations of ground filtering via PDAL, TIN-linear and Laplace interpolation via startin, 
natural neighbour (NN) interpolation via CGAL and radial IDW via GDAL/PDAL. I am currently working on integrating an ellipsoidal/quadrant-based IDW solution, but it is
proving much more difficult than anticipated (more about this at the end of the readme).

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
		* PDAL-IDW
    4. output format, one of:
        * ASC
        * GeoTIFF
    5. IDW interpolation radius in metres. _(optional)_
	6. IDW interpolation power in metres (the exponent of the inverse weighting). _(optional)_
	7. IDW interpolation fallback window size (more on this below). _(optional)_

An example call in the Windows Anaconda Prompt would be: `python C:/Users/geo-geek/some_folder/ip_main.py C:/Users/geo-geek/target_folder/ startin-TINlinear ASC`

Or for the IDW algorithm with radius and power values it would be `python C:/Users/geo-geek/some_folder/ip_main.py C:/Users/geo-geek/target_folder/ PDAL-IDW 10 2 GeoTIFF`

### A word of caution
If you are using an Anaconda virtual environment for PDAL/CGAL, you should first activate the environment in Anaconda prompt and _then_ run the relevant script
from the same prompt. So, for example:
1. `conda activate [environment_name]`
2. `python [file_path_to_main] [argument_1] [argument_2] [...]`

Relative file paths won't work in virtual environments, so make sure you specify the target folder using a full (absolute) file path.

**Note:** ASC export is not currently supported for the PDAL-IDW algorithm.

### Future work
The current version of this implementation runs as many processes in parallel as there are input files, hence using all processor cores when there are at least as many files as there are cores in
the given system. This multiprocessing implementation is based on Python built-in multiprocessing pools. A queue-based implementation would probably work better, but this is something for the scaling group to look at.

We should probably also experiment around with the `filters.pmf` method that is provided by PDAL. The example parametrisation uses `filters.smrf`. You can start experimenting with this simply by changing the JSON parametrisation, the code does not need to be edited.
There's further guidance in `gf_main.py` on what's what in `config.json`.

The PDAL-IDW workflow is actually built on top of GDAL, but since GDAL does not play well with Python data structures, I used the interface that is provided within PDAL's pipeline framework to implement it.
No part of the program currently uses the Python bindings of GDAL directly, but we might need to eventually start working with it. The ellipsoidal IDW features cannot be accessed through PDAL's interface for GDAL,
hence it cannot be used (hence PDAL-IDW only accepts one radius). There is a neat extra feature in the PDAL interface though, it allows a fallback method to be used. If you specify a value for an interpolation window
(the last command line argument above), wherever radial IDW fails, the program will look for values within a square kernel of pixels around the pixel that is being interpolated (presumably after the first round of
true IDW interpolation). For example, if you provide a value of 10 for this argument, it will look for values in a 10x10 square kernel around the pixel for values, weighting them based on their distance from the
pixel that is being interpolated. This can theoretically make the result more or less continuous (like the Voronoi and TIN-based methods).

I actually managed to implement an program that uses the GDAL interface directly (see the file `GDAL_attempt.py`) which is ready to read OGR vector files (for example ESRI shapefiles) and interpolate them
using GDAL's Python bindings. This would give access to the full IDW functionality of GDAL, but there is a problem: I implemented this because I though we could simply set the PDAL ground filtering implementation
to export into OGR files rather than LAS files and then simply use those as input for GDAL (GDAL is not compatible with LAS files). However, it turns out that there is a bug in the OGR writer of PDAL. It crashes
Python randomly while exporting. I can run the same code 10 times in a loop and 3-4 out of 10 attempts it will export correctly, but in the other cases it will crash. Unfortunately this makes it useless to us.
I have been trying to get this to work in Python 3.8, which could potentially be the source of the issues. I'll try to experiment with other combinations of Python and PDAL versions in the future.