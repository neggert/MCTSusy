from cmsscripts import das_utils, get_weights
import logging
import sys
from data_handling import *

print sys.argv

import json

input_json = sys.argv[1]
outfile=sys.argv[2]

with open(input_json) as f:
	params = json.load(f)

mctype = params['mctype']
mc_cat = params['mc_cat']
x_eff = params['x_eff']
files = map(str, params['files'])

logging.basicConfig(filename='/home/uscms33/'+mctype+".log", level=logging.INFO)

from itertools import *

def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

for group in grouper(10, files) :
	thesefiles = list(filter(None,group))
	logging.info( "Starting new file group" )
	logging.info( thesefiles[0])
	logging.debug( thesefiles )
	save_data_pandas( thesefiles, outfile, mctype, mc_cat, x_eff)
