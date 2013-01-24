#!/usr/bin/env python

#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    bkg_fit.py <signal_file> <mass1> <mass2> [-p] [--ncpu=<c>] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    --ncpu=<c>       Number of CPUs to use [default: 1]
    --channels=<c1,c2> Channels to use [default: of,sf]
    -p               Get p value (slow)

"""

from set_limits2 import *
from collections import defaultdict
import json

bkgs = ['of', 'vv', 'wjets', 'z']

def run_bonly_fit(sig_file, mass1, mass2, chans, ncpu, get_p):
    prefix = "limits/"+sig_file[:-5]+"_{}_{}_2".format(mass1, mass2)

    try:
        create_histfactory(sig_file, prefix, mass1, mass2, chans)
    except:
        return None

    rfile = R.TFile(prefix+"_combined_meas_model.root")

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    model.GetParametersOfInterest().first().setVal(0.)
    model.GetParametersOfInterest().first().setConstant()

    # run the fit
    res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save())

    fitPars = res.floatParsFinal()

    fitresults = {}
    for b in bkgs:
        fitvar = fitPars.find('n_{}'.format(b))
        fitresults[b] = (fitvar.getVal(), fitvar.getError())

    f = open("fit_results2.json", 'w')

    json.dump(fitresults, f, indent=3)

    f.close()


    model.SetSnapshot(model.GetParametersOfInterest())

    R.gROOT.ProcessLineSync(".L KS/AndersonDarlingTestStat.cc+")

    AD = R.RooStats.AndersonDarlingTestStat(model.GetPdf())

    # get the test statistic on data
    ts = AD.Evaluate(data, model.GetParametersOfInterest())

    if get_p:

        sampler = R.RooStats.ToyMCSampler(AD, 500)
        sampler.SetPdf(model.GetPdf())
        sampler.SetObservables(model.GetObservables())
        sampler.SetGlobalObservables(model.GetGlobalObservables())
        sampler.SetParametersForTestStat(model.GetParametersOfInterest())

        params = R.RooArgSet()
        params.add(model.GetNuisanceParameters())
        params.add(model.GetParametersOfInterest())

        if ncpu > 1:
            pc = R.RooStats.ProofConfig(ws, ncpu, "")
            sampler.SetProofConfig(pc)

        sampDist = sampler.GetSamplingDistribution(params)

        p = 1-sampDist.CDF(ts)

        print "P value:", p
        print "Test statistic on data: {:.7f}".format(ts)

        plot = R.RooStats.SamplingDistPlot()
        plot.AddSamplingDistribution(sampDist)

        plot.Draw()
        raw_input("...")

    print "Test statistic on data: {:.7f}".format(ts)

    return fitresults



if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)

    m1 = int(args['<mass1>'])
    m2 = int(args['<mass2>'])

    sig_file = args['<signal_file>']

    try:
        chans = args['--channels'].split(",")
    except AttributeError:
        chans = ['of', 'sf']

    ncpu = int(args['--ncpu'])

    get_p = bool(args['-p'])

    res = run_bonly_fit(sig_file, m1, m2, chans, ncpu, get_p)

