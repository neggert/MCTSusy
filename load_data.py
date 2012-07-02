from cmsscripts import das_utils
import sys
from data_handling import *

files = []

inputdslist = open(sys.argv[1], 'r')
datasets = inputdslist.read().splitlines()
inputdslist.close()

print datasets

files = []
prefix = "root://osg-se.cac.cornell.edu//xrootd/path/cms"
for ds in datasets :
    files.extend(das_utils.get_files_from_dataset(ds))

files = [prefix+f for f in files]

save_data_pandas( files, "/home/uscms33/data.hdf5", mctype="data")
