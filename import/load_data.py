from cmsscripts import das_utils, get_weights
import logging
import sys
from data_handling import *

print sys.argv

files = []

ds = sys.argv[1]


outfile=sys.argv[2]
mctype=sys.argv[3]
mc_cat = sys.argv[4]

try :
	xsec=float(sys.argv[5])
except IndexError :
	xsec = 1.

# get effective cross-section
x_eff = get_weights.get_eff_xsec( ds, xsec)

logging.basicConfig(filename='/home/uscms33/'+mctype+".log", level=logging.INFO)

logging.info( "Dataset: "+ds )
logging.info( "Output File"+outfile )
logging.info( "MCType"+mctype )

files = []
prefix = "root://osg-se.cac.cornell.edu//xrootd/path/cms"
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
	save_data_pandas( thesefiles, outfile, mctype, mc_cat, x_eff)
