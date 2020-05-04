# Ground Filtering, Interpolation, Hole Filling and Hydro-flattening

## In repo so far:
* `gf_processing.py`
* `gf_main.py`
* target_example (folder)
	* `config.json`
	* `fnames.txt`

The only working testing environment as of now is for pdal-based ground filtering via a rudimentary pool-based multiprocessing implementation.

## PDAL-based ground filtering (test version)
### User guide
1. Open `gf_main.py` and modify `target_folder` to point to the folder where the LAS files are.
2. Copy the example configuration files (`config.json`, `fnames.txt`) to your target folder.
3. Edit the configuration files.
	1. Specify the names of the LAS files in `fnames.txt` as shown in the example.
	2. Modify the configuration in `config.json`. Fine-tune the parametrisation to work well with the given data set.
4. Run `gf_main.py` **from the console**. If you run it from an IDE, it will probably not fork the processes properly.

If your pdal installation is in a virtual environment, the correct way to run this in `conda` is to issue these commands:
* `conda activate myenv`
* `python [file_path_to_main]`
Relative file paths won't work in virtual environments, so make sure you specify the target folder using a full (absolute) file path in `gf_main.py`.

### Future work
Keep in mind that this program can only fork as many processes as you have cores. On a hexacore processor, it will fork six processes and hence you are expected to specify six file names in `fnames.txt`.
It will ignore extra files. This could be improved via a queue-based framework, but this is merely a proof of concept and queues will probably be developed by the scaling group.

Furthermore, we should probably also experiment around with the `filters.pmf` method that is provided by GDAL. The example parametrisation uses `filters.smrf`. You can start experimenting with this simply by changing the JSON parametrisation, the code does not need to be edited.
There's further guidance on what's what in `config.json`, inside `gf_processing.py`.