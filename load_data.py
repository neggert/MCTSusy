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

from itertools import *

def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

for group in grouper(100, files) :
	thesefiles = list(filter(None,group))
	print thesefiles[0]
	save_data_pandas( thesefiles, "/home/uscms33/dataLoose.hdf5", mctype="data")
