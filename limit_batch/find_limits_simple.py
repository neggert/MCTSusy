#!/usr/bin/env python
"""
Combine hypothesis test results

Usage:
    combine_hypotest.py <output_file> <input_files>...

"""

import sys
import re
import ROOT as R
import docopt

args = docopt.docopt(__doc__)

files = args['<input_files>']

tests = {}

for f in files:
    # m1, m2 = map(int, re.search("_(\d+)_(\d+)_", f).groups())
    m1, m2 = map(int, re.search("_(\d+)_(\d+)[_\.]", f).groups())
    rf = R.TFile(f)
    if (m1,m2) not in tests.keys():
        tests[(m1,m2)] = rf.Get("result_sig_strength")
    else:
        tests[(m1,m2)].Add(rf.Get("result_sig_strength"))

with open(args['<output_file>'], 'w') as f:
    for masses in tests.keys():
        res = tests[masses]
        res.SetInterpolationOption(R.RooStats.HypoTestInverterResult.kSpline)
        m1,m2 = masses
        exp = res.GetExpectedUpperLimit()
        exp_down = res.GetExpectedUpperLimit(-1)
        exp_up = res.GetExpectedUpperLimit(+1)
        obs = res.UpperLimit()

        f.write("{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4:.2f}\t{5:.2f}\n".format(m1, m2, exp, exp_down, exp_up, obs))



