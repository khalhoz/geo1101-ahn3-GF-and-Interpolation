# Ground Filtering, Interpolation, Hole Filling and Hydro-flattening Testing Environment

**This is the Python repo of the ground filtering/interpolation team of the AHN3 GEO1101 (synthesis project) group of 2020.**

## In repo so far:

* `Classification_F1_score` _(Khaled & Manos's quantitative ground filtering evaluation results)_
* `Khaled_Interpolation` _(Some interpolation testing code from Khaled)_
* `gf_target_example` _(example folder for LAS input files)_
	* `config.json` _(ground filtering configuration file which you can push updates to)_
	* `config_default.json` _(default ground filtering configuration file which should not be modified)_
	* `config_preprocess.json` _(Khaled & Manos's final pipeline - only pre-processes, does not ground-filter)_
	* `fnames.txt` _(list of input files names for the ground filtering program)_
* `Ground filtering.docx` _(Manos's draft report on the ground filtering results)_
* `README.md` _(this readme file)_
* `gdal_attempt.py` _(a half-baked test program for direct interpolation via GDAL)_
* `gf_main.py` _(main file for ground filtering)_
* `gf_processing.py` _(ground filtering code)_
* `ip_main.py` _(main file for interpolation)_
* `ip_processing.py` _(interpolation code)_
* `ip_st_test` _(some raster tile statistics generation code from Khaled)_
* `las_prepare.py` _(factored out from `ip_processing.py`)_
* `split_bbg_into_subfiles.py` _(Lisa's code for splitting the [BGG](https://www.pdok.nl/introductie/-/article/cbs-bestand-bodemgebruik) into subsets of river, sea, and other water body polygons)_

The testing environment so far includes multiprocessing pool-based implementations of ground filtering/pre-processing via PDAL, TIN-linear and Laplace interpolation via startin, 
constrained Delaunay-based (CDT) TIN-linear and natural neighbour (NN) interpolation via CGAL, radial IDW via GDAL/PDAL and quadrant-based IDW via scipy cKDTree and our own code.
We might be able to implement an ellipsodial IDW solution in the future, but it is proving much more difficult than anticipated (more about this at the end of the readme).

**NEW STUFF**

* First release of CGAL-CDT. Everything seems to be working, but needs more testing.
* Added an option to pre-process the data before interpolation using a PDAL pipeline, creating a single entry point for the whole process
* Re-wrote the CGAL-NN algorithm to use true natural neighbours and not regular neighbours
* Fixed lots of bugs

## Current tasks

I am working with Lisa on hole-filling and hydro-flattening at the moment. For the final production environment, we aim to skip secondary hole filling altogether by using a TIN-based
interpolation method. This will lend us some extra time to develop a useful water body and river flattening algorithm. We already have the necessary shapefiles from BBG thanks to Lisa, and
will now move on to implementing algorithms. Optimally, we are hoping to make the river polygons part of the triangulation via the CDT implementation and use the hole boundaries for water body
shore elevation values (which will be used to set the interpolated value inside the holes that arise due to the constraints). As a backup plan, we would use Laplace interpolation, but in this
case we would not be able to create holes in the triangulation using the water polygons, and we would probably end up needing to overlay them with the rasters post-interpolation.

Following some further consultation with Jeroen, we will also double-check how well our primary interpolation method interpolates the building gaps. If they look too triangular, we might
consider trying to flatten them (like water bodies), or using an alternative secondary interpolation mechanism inside of them.

## PDAL pipeline entry point (stand-alone ground filtering/pre-processing)

This is a secondary entry point that only runs a PDAL pipeline. It is intended for testing PDAL's ground filtering and pre-processing capabilities.

1. Copy the example configuration files (`config.json`, `fnames.txt`) to your target folder (in which the LAS files are located).
2. Edit the configuration files.
	1. Specify the names of the LAS files in `fnames.txt` as shown in the example.
	2. Modify the pdal pipeline configuration in `config.json`. Fine-tune the parametrisation to work well with the given data set.
4. Run `gf_main.py` **from the console**. If you run it from an IDE, it will probably not fork the processes properly. It takes one command line argument,
which should be the target folder. An example call in the Windows Anaconda Prompt would be: `python C:/Users/geo-geek/some_folder/gf_main.py C:/Users/geo-geek/target_folder/`

It will deposit the ground filtered tiles as LAS files tagged with `_gf.las` in the target folder.
You can provide an additional CMD argument to modify the default output tag (which is `_gf`). You don't need the underscore, to export `somefile_dsm.las` for example, you may write
`python [file_path_to_main] [target_folder] dsm`. You can also provide a tag for the input config file. For example for a pre-processing job you may want to call it `config_pre.json`,
in which case you would use `python [file_path_to_main] [target_folder] dsm pre`.

## Primary entry point (optional ground filtering/pre-processing + interpolation)

`ip_main.py` is intended to be the primary entry point to the program. You can either use it to interpolate pre-processed LAS files you had already generated using the above secondary entry point,
or run ground-filtering/pre-processing as part of this process. In the prior case, the program assumes that you used the default `_gf.las` tag for your final exports. It also needs the same
`fnames.txt` file to be present (and the file names in it shoud **not** be tagged `_gf.las`).

If you use this as the primary entry point, you do not need to separately pre-process the data. Just provide `True` as the second argument (as described below), and the program will automatically
run the PDAL pipeline too. In this case, it is necessary that your config JSON is called `config_preprocess.json`. Note that this solution transfers data between the PDAL-based preprocessing and
the interpolation in memory, and is hence _faster_ than running `gf_main.py` separately beforehand. The PDAL-IDW method has also been adapted to work so that pre-processing an interpolation form
a single pipeline if you are also running pre-processing from this entry point.

The intended workflow is:
1. Check that you still have the ground filtered files and the `fnames.txt` file in your target folder if you had run the ground filtering/pre-processing workflow separately beforehand.
2. Run `ip_main.py` **from the console**. If you run it from an IDE, it will probably not fork the processes properly. The following arguments should be provided:
    1. target folder (most likely the same as the one you used with PDAL)
	2. bool to indicate if you would like pre-processing to be run before interpolation _(the default value is `False`)_
    3. pixel size (in metres) for interpolation _(the default value is 1)_
    4. interpolation method, one of:
        * startin-TINlinear
        * startin-Laplace _(default)_
        * CGAL-NN _(NOTE: Re-designed the algorithm, now uses true natural neighbours.)_
		* CGAL-CDT
		* PDAL-IDW
		* IDWquad
	5. output format, one of:
        * ASC
        * GeoTIFF _(default)_
    6. IDW argument 1 / File path for input polygons for CDT:
		* _If using PDAL-IDW:_ IDW interpolation radius in metres
		* _If using IDWquad:_ The _starting_ radius/number of neighbours _k_ to query
		* _If using CGAL-CDT:_ A file path to a fiona-compatible vector file to be used for setting up the CDT constraints.
	7. IDW argument 2: IDW interpolation power (exponent) in metres _(both for PDAL-IDW and IDWquad)_
	8. IDW argument 3:
		* _If using PDAL-IDW:_ interpolation fallback window size
		* _If using IDWquad:_ minimum number of points to find per quadrant
	9. IDW argument 4: query radius/number of neighbours _k_ to query, increment step value _(only for IDWquad)_
	10. IDW argument 5: IDWquad method, one of:
		* radial _(for iterative radius increments)_
		* k-nearest _(for iterative increments of how many neighbours to query)_
	11. IDW argument 6: IDWquad KD-tree query tolerance value _eps_
	12. IDW argument 7: IDWquad maximum number of iterations before declaring no-data and proceeding to next pixel

An example call in the Windows Anaconda Prompt would be:

`python C:/Users/geo-geek/some_folder/ip_main.py C:/Users/geo-geek/target_folder/ False 2 startin-TINlinear ASC`

Or for the PDAL-IDW algorithm (using pre-processing) with radius and power values it would be

`python C:/Users/geo-geek/some_folder/ip_main.py C:/Users/geo-geek/target_folder/ True 0.5 PDAL-IDW GeoTIFF 10 2`

As indicated above, for CGAL-CDT the sixth parameter is a file path too. For example:

`python C:/Users/geo-geek/some_folder/ip_main.py C:/Users/geo-geek/target_folder/ False 1.5 CGAL-CDT GeoTIFF C:/Users/geo-geek/some_folder/killer_polygons.shp`

## A word of caution

If you are using an Anaconda virtual environment for PDAL/CGAL, you should first activate the environment in Anaconda prompt and _then_ run the relevant script
from the same prompt. So, for example:
1. `conda activate [environment_name]`
2. `python [file_path_to_main] [argument_1] [argument_2] [...]`

Relative file paths won't work in virtual environments, so make sure you specify the target folder using a full (absolute) file path.
I implemented the program in such a way, that only those modules are imported which are necessary for the given operation. For instance, if you are using neither pre-processing,
nor PDAL-IDW, then you do not need to have PDAL in your Python environment to run the program.

Another word of caution with the outputs is that they all use a fixed no-data value of -9999. This includes the GeoTIFF exporter. To view the results correctly, you should keep in
mind that while the upper bounds of the data will be determined correctly by the viewer software (e.g. QGIS), the lower bound will be -9999. To see the DTM/DSM the program interpolated,
you need to set the lower bound of the colour scale to a higher value relevant to the data. As AHN3 depicts Dutch terrain, you are advised to use a value somewhere between -20 and 0 metres
depending on the tile. Negative elevation values are perfectly possible in The Netherlands.
In QGIS, you do this by right clicking on your raster layer, and clicking on "Properties...". In the window that pops up, you can change the lower bound of the colour scale by
adjusting the value in the field Symbology --> Band Rendering --> Min.

**Note:** ASC export is not currently supported for the PDAL-IDW algorithm.

**Another note:** You are advised to configure the IDWradial parametrisation **with performance in mind** when first getting started with IDWquad. Otherwise it might take _veeeeeery long_ to finish.

## More about the IDW algorithms

### PDAL-IDW
The PDAL-IDW workflow is actually built on top of GDAL, but since GDAL does not play well with Python data structures, I used the interface that is provided within PDAL's pipeline framework to
implement it. No part of the program currently uses the Python bindings of GDAL directly, but we might need to eventually start working with it. The ellipsoidal IDW features cannot be accessed
through PDAL's interface for GDAL, hence they cannot be used here (hence PDAL-IDW only accepts one radius). There is a neat extra feature in the PDAL interface though, it allows a fallback method
to be used. If you specify a value for an interpolation window (IDW argument 3 above), wherever radial IDW fails, the program will look for values within a square kernel of pixels around the pixel
that is being interpolated (presumably after the first round of true IDW interpolation). For example, if you provide a value of 10 for this argument, it will look for values in a 10x10 square
kernel around the pixel for values, weighting them based on their distance from the pixel that is being interpolated. This can theoretically make the result more or less continuous (a bit more
like the Voronoi and TIN-based methods).

### IDWquad
This is a quadrant-based IDW implementation that is not built on top of third-party software (apart from scipy, from which cKDTree is used). It builds a KD-tree representation of the points of the
input tile and overlay it with a raster of the desired dimensions. For each pixel, it iteratively queries more and more points until it has enough points **per quadrant** to consider an IDW
interpolation reliable. The algorithm can either based its KD-tree queries on k-nearest neighbours to find, or a query radius to search within.
I'll explain how it works by giving some more detail about the parameters, listing them in the same order as in the list of arguments above.

* **Starting interpolation radius/number of neighbours _k_ to query:** The starting value for the radius or number of neighbours to query. This is the value that will be incremented until enough points per quadrant can be found to interpolate a value for the pixel.
* **Interpolation power:** this is self-explanatory, it is the power that is used _if_ IDW interpolation can take place based on the number of neighbours per quadrant that were found.
* **Minimum number of points to find per quadrant:** the minimum number of neighbours to find in _each_ quadrant to interpolate a value for the pixel. E.g. if it is set to one, the algorithm will be satisfied by finding only one neighbour per quadrant. If you set it to a higher value, it will increase the radius/number of neighbours to find until it finds that many neighbours per quadrant (or reaches the maximum number of iterations).
* **Increment value:** the algorithm will increase the radius/number of closest neighbours queried by this amount in each iteration until there are enough neighbours per quadrant to interpolate.
* **Method:** as already explained, you can either specify "radial" to make the algorithm iteratively expand its search radius, or "k-nearest" to iteratively increase the number of nearest neighbours to query.
* **Tolerance value _eps_:** this is very important in terms of performance. As we discussed with Ravi on Discord, the best way to increase performance here is to use _approximate_ KD-tree queries rather than exact queries. Fine tune this parameter to drastically improve performance. This does not affect continuity, it does however affect quality. You can find more info about this on the following pages:
	* If you use radial queries: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_ball_point.html#scipy.spatial.cKDTree.query_ball_point
	* If you use k-nearest queries: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html#scipy.spatial.cKDTree.query
* **Iteration limit:** you may further fine-tune performance by setting a limit on how many times the algorithm may increment the radius/number of neighbours to find. In some cases this may also drastically improve performance at the cost of continuity.

## Future work

The current version of this implementation runs as many processes in parallel as there are input files, hence using all processor cores when there are at least as many files as there are cores in
the given system. This multiprocessing implementation is based on Python built-in multiprocessing pools. A queue-based implementation would probably work better, but this is something for the scaling group to look at. _No final decision yet._

We should probably also experiment around with the `filters.pmf` method that is provided by PDAL. The example parametrisation uses `filters.smrf`. You can start experimenting with this simply by changing the JSON parametrisation, the code does not need to be edited.
There's further guidance in `gf_main.py` on what's what in `config.json`. _(Completed.)_

I actually managed to implement an program that uses the GDAL interface directly (see the file `GDAL_attempt.py`) which is ready to read OGR vector files (for example ESRI shapefiles) and interpolate them
using GDAL's Python bindings. This would give access to the full IDW functionality of GDAL, but there is a problem: I implemented this because I though we could simply set the PDAL ground filtering implementation
to export into OGR files rather than LAS files and then simply use those as input for GDAL (GDAL is not compatible with LAS files). However, it turns out that there is a bug in the OGR writer of PDAL. It crashes
Python randomly while exporting. I can run the same code 10 times in a loop and 3-4 out of 10 attempts it will export correctly, but in the other cases it will crash. Unfortunately this makes it useless to us.
I have been trying to get this to work in Python 3.8, which could potentially be the source of the issues. I'll try to experiment with other combinations of Python and PDAL versions in the future.
_Now that I also implemented quadrant-based IDW, perhaps we might not need this after all?_

A further idea would be to implement no-data mask layers for the GeoTIFF exporter. This would make them easier to visualise, as the no-data pixels would be masked out by the no-data
layer, rather than simply being marked as no-data by the fixed value -9999. _Not being worked on at the moment._