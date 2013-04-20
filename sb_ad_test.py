#!/usr/bin/env python

#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    bkg_fit.py <filename> [--ncpu=<c>] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    --ncpu=<c>       Number of CPUs to use [default: 1]
    --channels=<c1,c2> Channels to use [default: of,sf]

"""

from set_limits import *
from collections import defaultdict
import json
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import IndexLocator, FixedFormatter
import numpy as np

bkgs = ['top', 'vv', 'wjets', 'z']

def run_bonly_fit(file_name, ncpu, data_prefix="data", data_file_name="data.root"):

    rfile = R.TFile(file_name)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    sbmodel = ws.obj("ModelConfig")
    bmodel = sbmodel.Clone()
    bmodel.SetName("bkgModelConfig")

    constr = sbmodel.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    bmodel.GetParametersOfInterest().first().setVal(0.)
    bmodel.GetParametersOfInterest().first().setConstant()
    bmodel.SetSnapshot(bmodel.GetParametersOfInterest())    


    pars = sbmodel.GetNuisanceParameters()

    # run the fit
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR)
    bmodel.LoadSnapshot()
    res = bmodel.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.PrintLevel(0))

    sbmodel.GetParametersOfInterest().first().setVal(1.)
    sbmodel.GetParametersOfInterest().first().setConstant()
    # sbmodel.SetSnapshot(sbmodel.GetParametersOfInterest())    

    params = R.RooArgSet()
    params.add(sbmodel.GetNuisanceParameters())
    params.add(sbmodel.GetParametersOfInterest())
    sbmodel.SetSnapshot(params)

    # getattr(ws, "import")(sbmodel)
    getattr(ws, "import")(bmodel)


    results = []
    from IPython.parallel import Client

    rc = Client()
    dview = rc[:]
    dview.execute("import ROOT as R")
    dview.execute("R.gROOT.ProcessLineSync('.L KS/AndersonDarlingTestStat.cc+')")
    lview = rc.load_balanced_view()

    for i in xrange(50):
        r = lview.apply_async(get_sig_p_value, ws, 10)
        results.append(r)

    lview.wait(results)

    sigSampleDist = R.RooStats.SamplingDistribution()
    for r in results:
        sigSampleDist.Add(r.result)


    f = R.TFile("BkgADDist.root")
    bkgSampleDist = f.Get("sampDist")

    c1 = R.TCanvas()
    plot = R.RooStats.SamplingDistPlot(50, 0, 10)
    plot.AddSamplingDistribution(sigSampleDist)
    plot.AddSamplingDistribution(bkgSampleDist)

    plot.Draw()
    c1.SaveAs("plots/adtest.pdf")


    # raw_input("...")


def get_sig_p_value(ws, n):
    sbmodel = ws.obj("ModelConfig")
    bmodel = ws.obj("bkgModelConfig")

    R.gROOT.ProcessLineSync(".L KS/AndersonDarlingTestStat.cc+")

    AD = R.RooStats.AndersonDarlingTestStat(bmodel.GetPdf())

    sampler = R.RooStats.ToyMCSampler(AD, n)
    sampler.SetPdf(sbmodel.GetPdf())
    sampler.SetObservables(sbmodel.GetObservables())
    sampler.SetGlobalObservables(sbmodel.GetGlobalObservables())
    bparams = bmodel.GetSnapshot()
    sampler.SetParametersForTestStat(bmodel.GetSnapshot())
    sampler.SetSamplingDistName("Signal+Background")

    sigSampleDist = sampler.GetSamplingDistribution(sbmodel.GetSnapshot())

    return sigSampleDist


if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)

    file_name = args['<filename>']



    ncpu = int(args['--ncpu'])


    res = run_bonly_fit(file_name, ncpu)

