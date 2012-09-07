from cmsscripts import das_utils
import logging
import sys
from data_handling import *


files = []

inputdslist = open(sys.argv[1], 'r')
datasets = inputdslist.read().splitlines()
inputdslist.close()


outfile=sys.argv[2]
mctype=sys.argv[3]
try :
	weight=sys.argv[4]
except IndexError :
	weight = 1.

logging.basicConfig(filename='/home/uscms33/'+mctype+".log", level=logging.INFO)

logging.info( datasets )
logging.info( outfile )
logging.info( mctype )

files = []
prefix = "root://osg-se.cac.cornell.edu//xrootd/path/cms"
for ds in datasets :
    files.extend(das_utils.get_files_from_dataset(ds))

logging.info( str(len(files))+ " files")

files = [prefix+f for f in files]

logging.debug(files)

from itertools import *

def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

for group in grouper(100, files) :
	thesefiles = list(filter(None,group))
	logging.info( "Starting new file group" )
	logging.info( thesefiles[0])
	logging.debug( thesefiles )
	save_data_pandas( thesefiles, outfile, mctype=mctype, weight=weight)
