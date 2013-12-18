#!/usr/bin/env python
"""
Combine hypothesis test results

Usage:
    combine_hypotests.py <output_prefix> <input_files>... 

"""

import sys
import re
import ROOT as R
import docopt
import scipy.interpolate
import scipy.optimize
import scipy.stats
import sklearn.gaussian_process as gp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from copy import deepcopy
import os.path

import warnings
warnings.simplefilter("ignore", DeprecationWarning)



def main(input_files, output_prefix):

    tests = {}

    for f in input_files:
        m1, m2 = map(int, re.search("_(\d+)_(\d+)_", f).groups())
        rf = R.TFile(f)
        if (m1,m2) not in tests.keys():
            tests[(m1,m2)] = deepcopy(rf.Get("result_sig_strength"))
        else:
            tests[(m1,m2)].Add(rf.Get("result_sig_strength"))
        rf.Close()

    for masses in tests.keys():
        res = tests[masses]
        m1,m2 = masses

        fname = output_prefix+"{0}_{1}.root".format(m1, m2)
        if os.path.exists(fname):
            rf = R.TFile(fname)
            res.Add(rf.Get("result_sig_strength"))
            rf.Close()

        f = R.TFile(output_prefix+"{0}_{1}.root".format(m1, m2), "RECREATE")
        f.cd()
        res.Write()
        f.Close()        

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    files = args['<input_files>']

    output_prefix = args['<output_prefix>']
    main(files, output_prefix)

