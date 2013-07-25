#!/usr/bin/env python

#! /usr/bin/env python

"""Evaluate p-value for a specific model

Usage:
    sensitivity.py <signal_file> <mass1> <mass2> [-h][--channels=<c1,c2>]
    sensitivity.py batch <signal_file> <mass_file> <jobnum> <output_file> [-h] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    --channels=<c1,c2> Channels to use [default: of,sf]

"""

from set_limits import *
from collections import defaultdict
import json

bkgs = ['top', 'vv', 'wjets', 'z']

def frequentist_pvalue(filename):
    rfile = R.TFile(filename)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    sbmodel = ws.obj("ModelConfig")

    sbmodel.GetParametersOfInterest().first().setVal(1.)
    sbmodel.GetParametersOfInterest().first().setConstant(True)

    constr = sbmodel.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    # run the initial fit
    sbmodel.GetPdf().fitTo(data, R.RooFit.Constrain(constr))

    poi_hat = sbmodel.GetParametersOfInterest().first().getVal()
    poi_hat_err = sbmodel.GetParametersOfInterest().first().getError()

    sbmodel.SetSnapshot(sbmodel.GetParametersOfInterest())

    bmodel = sbmodel.Clone()

    bmodel.GetParametersOfInterest().first().setVal(0.)
    bmodel.GetParametersOfInterest().first().setConstant(True)
    bmodel.SetName("ModelConfig_bonly")
    bmodel.SetSnapshot(bmodel.GetParametersOfInterest())

    calc = R.RooStats.FrequentistCalculator(data, sbmodel, bmodel)
    calc.SetToys(200, 100)

    R.gROOT.ProcessLineSync(".L KS/AndersonDarlingTestStat.cc+")
    AD = R.RooStats.AndersonDarlingTestStat(bmodel.GetPdf())

    # get the test statistic on data
    # ts = AD.Evaluate(data, model.GetParametersOfInterest())

    sampler = calc.GetTestStatSampler()
    sampler.SetNEventsPerToy(data.numEntries())
    sampler.SetGenerateBinned(True)
    sampler.SetTestStatistic(AD)

    htr = calc.GetHypoTest()
    htr.SetPValueIsRightTail(True)
    htr.SetBackgroundAsAlt(False)
    htr.Print()

def run_limit(sig_file, mass1, mass2, chans):
    prefix = "limits/"+sig_file[:-5]+"_{0}_{1}".format(mass1, mass2)

    # try:
    create_histfactory(sig_file, prefix, mass1, mass2, chans)
    # except:
        # return None

    frequentist_pvalue(prefix+"_combined_meas_model.root")

if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)

    if not args['batch']:
        m1 = int(args['<mass1>'])
        m2 = int(args['<mass2>'])
    else:
        with open(args['<mass_file>']) as f:
            masses = json.load(f)
        m1, m2 = [int(m) for m in masses[int(args['<jobnum>'])]]

    sig_file = args['<signal_file>']

    try:
        chans = args['--channels'].split(",")
    except AttributeError:
        chans = ['of', 'sf']

    run_limit(sig_file, m1, m2, chans)


