### MAIN FILE FOR PDAL TESTING ###

from gf_processing import initialise, start_pool

# You need to specify the folder in which you LAS files are.
# The folder needs to contain the config JSON in the right format (see example file), with name "config.json".
# The folder also needs to contain a file called "fnames.txt" with the names of the files you want to process,
# as also shown in the example file.
# Use an ABSOLUTE file path to avoid problems (especially when working from virtual environments)!
target_folder = "C:/Users/geo-geek/example_folder/"
# The output files will be written to the target folder too. They will be tagged with "_out".

def main():
    config, fnames = initialise(target_folder)
    start_pool(target_folder, config, fnames)
    print("Success.")

if __name__ == '__main__':
    main()