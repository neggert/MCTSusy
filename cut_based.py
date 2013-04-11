#!/usr/bin/env python

#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    cut_based.py <filename> 

Options:
    -h --help        Show this screen.

"""

# from set_limits import *
from collections import defaultdict
import json
import matplotlib.pyplot as plt
import numpy as np
# import ROOT as R
# plt.switch_backend("pdf")
from copy import deepcopy

from IPython.parallel import Client

rc = Client()
dview = rc[:]
with dview.sync_imports():
    import cut_based_libs
    import ROOT

def run_cut_based(file_name, ncpu):

    rfile = ROOT.TFile(file_name)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    ROOT.RooStats.RemoveConstantParameters(constr)

    model.GetParametersOfInterest().first().setVal(0.)
    model.GetParametersOfInterest().first().setConstant()

    obs_sf = ws.obj("obs_x_sf")
    obs_of = ws.obj("obs_x_of")

    obs_sf.setRange("fitRange", 10., 120.)
    obs_of.setRange("fitRange", 10., 120.)


    # run the fit
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, ROOT.RooFit.Constrain(constr), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(0),
                               ROOT.RooFit.Range("fitRange"))

    params = model.GetParametersOfInterest()
    params.add(model.GetNuisanceParameters())
    model.SetSnapshot(params)


    data_results = cut_based_libs.get_results(ws, res)
    print json.dumps(data_results, indent=4)

    n_events = model.GetPdf().getPdf("of").expectedEvents(ROOT.RooArgSet(obs_of))+\
                model.GetPdf().getPdf("sf").expectedEvents(ROOT.RooArgSet(obs_sf))

    n_toy_events = 10000
    true = data_results['sf']['sum']['high'][0]*n_toy_events/n_events
    # local test
    # print fit_toy(ws, 1000000)

    lview = rc.load_balanced_view()

    means = []
    widths=[]

    n = [10000, 50000, 70000, 100000, 200000, 300000, 400000, 500000, 700000, 1000000]

    all_results = defaultdict(list)

    for n_toy_events in n:
        results = []
        true = data_results['sf']['sum']['high'][0]*n_toy_events/n_events
        # true = data_results['sf']['sum']['high'][0]


        for i in xrange(10):
            r = lview.apply_async(fit_toy, deepcopy(ws), n_toy_events)
            results.append(r)

        lview.wait(results)

        sf_vals = []
        sf_errs = []
        for r in results:
            all_results[n_toy_events].append(r.result)
            sf_vals.append(r.result['sf']['sum']['high'][0])
            sf_errs.append(r.result['sf']['sum']['high'][1])

        sf_vals = np.asarray(sf_vals)
        sf_errs = np.asarray(sf_errs)

        pulls = (sf_vals-true)/sf_vals.std()

        means.append(pulls.mean())
        widths.append(pulls.std())

    means = np.asarray(means)
    widths = np.asarray(widths)
    plt.plot(n, means, color="k")
    # plt.fill_between(n, means-widths, means+widths, color="b", alpha=0.5)
    plt.show()


    import IPython
    IPython.embed()

def fit_toy(ws, n):
    model = ws.obj("ModelConfig")

    params = model.GetParametersOfInterest()
    params.add(model.GetNuisanceParameters())

    dummy = ROOT.RooStats.NumEventsTestStat(model.GetPdf()) # this doesn't do anything, we just need something to give ToyMCSampler
    mc = ROOT.RooStats.ToyMCSampler(dummy, 1)
    mc.SetPdf(model.GetPdf())
    mc.SetObservables(model.GetObservables())
    mc.SetGlobalObservables(model.GetGlobalObservables())
    mc.SetParametersForTestStat(model.GetParametersOfInterest())
    mc.SetNEventsPerToy(n)

    constr = model.GetNuisanceParameters()
    ROOT.RooStats.RemoveConstantParameters(constr)

    toy_data = mc.GenerateToyData(model.GetSnapshot())
    res = model.GetPdf().fitTo(toy_data, ROOT.RooFit.Constrain(constr), ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
    
    return cut_based_libs.get_results(ws, res)


if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)
    print args

    file_name = args['<filename>']



    # ncpu = int(args['--ncpu'])

    res = run_cut_based(file_name, 0)

