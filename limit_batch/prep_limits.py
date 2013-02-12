#!/usr/bin/env python
"""Prepare file used to steer limit setting

Usage:
    prep_limits.py <output_file> <npoints> <file>...

Arguments:
    <output_file>
    <npoints>   Number of points in poi space to scan for each model
    <file>...   List of input model files
"""

import docopt
import ROOT as R
import numpy as np

args = docopt.docopt(__doc__)

npoints = int(args['<npoints>'])
files = args['<file>']
output = args['<output_file>']

for filename in files:
    rfile = R.TFile(filename)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    sbmodel = ws.obj("ModelConfig")

    constr = sbmodel.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    # run the initial fit
    sbmodel.GetPdf().fitTo(data, R.RooFit.Constrain(constr))

    poi_hat = sbmodel.GetParametersOfInterest().first().getVal()
    poi_hat_err = sbmodel.GetParametersOfInterest().first().getError()

    poi_vals = np.linspace(max([0.001, poi_hat-2*poi_hat_err]), poi_hat+8*poi_hat_err, npoints)

    with open(output, "a") as out:
        for p in poi_vals:
            line = '\t'.join([filename, str(p)])
            out.write(line+"\n")
