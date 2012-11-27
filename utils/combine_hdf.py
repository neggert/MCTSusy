#!/usr/bin/env python

"""
This script combines the dataframe stored in the "data" key from all input files
into a single dataframe in a single output file.

Useage: ./combin_hdf.py input_glob output_file

Arguments:
    input_glob -- linux shell style pattern that can be matched to input files
    output_file -- ouptut hdf file
"""
import sys
from pandas import *
import os

if len(sys.argv) < 4:
    sys.exit("Need at least two input files and an output file")

output_file = sys.argv[-1]
if os.path.exists(output_file) :
    sys.exit("{} already exists. Exiting.".format(output_file))

stores = []

for f in sys.argv[1:-1]:
    print "Opening file {}".format(f)
    stores.append(HDFStore(f))

df = stores[0]['data']

for s in stores[1:]:
    df = df.append(s['data'], ignore_index=True)

df = df.drop_duplicates(['run', 'lumi', 'event'])

out = HDFStore(sys.argv[-1])
print "Saving to file {}".format(sys.argv[-1])
out.put("data", df)
out.flush()
out.close()


