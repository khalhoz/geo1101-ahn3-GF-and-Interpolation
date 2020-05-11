# Ground Filtering, Interpolation, Hole Filling and Hydro-flattening

## In repo so far:
* `gf_processing.py`
* `gf_main.py`
* target_example (folder)
	* `fnames.txt`
	* `config.json` _(which you can push updates to)_
	* `config_default.json` _(which should not be modified, it is for reference)_

The only working testing environment as of now is for pdal-based ground filtering via a rudimentary pool-based multiprocessing implementation.

## PDAL-based ground filtering (test version)
### User guide
1. Copy the example configuration files (`config.json`, `fnames.txt`) to your target folder (in which the LAS files are located).
2. Edit the configuration files.
	1. Specify the names of the LAS files in `fnames.txt` as shown in the example.
	2. Modify the pdal pipeline configuration in `config.json`. Fine-tune the parametrisation to work well with the given data set.
4. Run `gf_main.py` **from the console**. If you run it from an IDE, it will probably not fork the processes properly. It takes one command line argument, which should be the target folder. An example call in the Windows Anaconda Prompt would be: python `C:/Users/geo-geek/some_folder/gf_main.py C:/Users/geo-geek/target_folder/`

If your pdal installation is in a conda virtual environment, you need to first activate your virtual environment and _then_ run the program:
1. `conda activate [environment_name]`
2. `python [file_path_to_main] [target_folder]`

Relative file paths won't work in virtual environments, so make sure you specify the target folder using a full (absolute) file path.

### Future work
The current version of this implementation runs as many processes in parallel, as there are input files, hence using all processor cores. This multiprocessing implementation is based
on Python built-in multiprocessing pools. A queue-based implementation would probably work better, but this is something for the scaling group to look at.

Furthermore, we should probably also experiment around with the `filters.pmf` method that is provided by GDAL. The example parametrisation uses `filters.smrf`. You can start experimenting with this simply by changing the JSON parametrisation, the code does not need to be edited.
There's further guidance in `gf_main.py` on what's what in `config.json`.