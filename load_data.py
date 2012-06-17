from cmsscripts import das_utils
import sys
from data_handling import *

files = []

inputdslist = open(sys.argv[1], 'r')
datasets = inputdslist.read().splitlines()
inputdslist.close()

print datasets

files = []
for ds in datasets :
    files.extend(das_utils.get_files_from_dataset(ds))

save_data_pandas( files, "data.hdf5", mctype="data")
