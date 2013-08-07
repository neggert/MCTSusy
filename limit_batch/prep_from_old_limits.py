#!/usr/bin/env python
"""Prepare file used to steer limit setting

Usage:
    prep_from_old_limits.py <output_file> <npoints> <file>

Arguments:
    <output_file>
    <npoints>   Number of points in poi space to scan for each model
    <file>      Old limits to use as a starting point
"""

import docopt
import numpy as np

args = docopt.docopt(__doc__)

npoints = int(args['<npoints>'])
in_file = args['<file>']
out_file = args['<output_file>']

in_data = np.loadtxt(in_file)

with open(out_file, 'w') as f:
    for line in in_data:
        m1, m2, obs, exp, expm1, expp1 = line
        m1, m2 = int(m1), int(m2)
        lims = [obs, exp, expm1, expp1]
        low = min(lims)
        high = max(lims)
        center = (low+high)/2
        width = max([(high-low)/2, 0.5])
        bottom = max([center-3*width, 0.01])
        top = center+3*width

        vals = np.linspace(bottom, top, npoints)
        for v in vals:
            out_line = "../limits/chi_templates_{0}_{1}_combined_meas_model.root {2}\n".format(m1, m2, v)
            f.write(out_line)
