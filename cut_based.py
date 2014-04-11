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
import scipy.stats
# import ROOT as R
# plt.switch_backend("pdf")
from copy import deepcopy

from IPython.parallel import Client

import bkg_fit
import ROOT

rc = Client()
dview = rc[:]
with dview.sync_imports():
    import bkg_fit
    import ROOT
dview.execute("ROOT.gROOT.ProcessLineSync('.L IntegralErrorShape.C+')")

ROOT.gROOT.ProcessLineSync('.L IntegralErrorShape.C+')

def run_cut_based(file_name, ncpu):

    high = 200.
    rfile = ROOT.TFile(file_name)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    ROOT.RooStats.RemoveConstantParameters(constr)

    pars = model.GetNuisanceParameters()


    # run the fit
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    counting = False
    ws.obj("obs_x_sf").setRange("fit", 10, 120.)
    ws.obj("obs_x_of").setRange("fit", 10, 120.)
    res = model.GetPdf().fitTo(data, ROOT.RooFit.Constrain(constr), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Range("fit"))

    r = ROOT.get_results(ws, res, 120., high)
    data_results = bkg_fit.results_to_dict_counting(r)
    model.SetSnapshot(pars)


    # test
    print fit_toy(ws, 2)

    lview = rc.load_balanced_view()

    n_tries = 1
    all_results = [[] for _ in xrange(n_tries)]

    true_sum_highs_sf = []
    true_sum_highs_of = []





    for j in xrange(n_tries):

        model.LoadSnapshot()

        true_sum_highs_sf.append(data_results['sf']['sum']['high'][0])
        true_sum_highs_of.append(data_results['of']['sum']['high'][0])

        results = []
        
        for i in xrange(8):
            r = lview.apply_async(fit_toy, ws, 2)
            results.append(r)

        lview.wait(results)

        for r in results:
            all_results[j].extend(r.result)


        dview.results.clear()
        # dview.clear(block=True)

    import IPython
    IPython.embed()



    means_sf = np.asarray([np.median([x['sf']['sum']['high'][0] for x in a]) for a in all_results])
    std_sf = np.asarray([np.std([x['sf']['sum']['high'][0] for x in a]) for a in all_results])

    pull_means_sf = np.asarray([np.mean( [(x['sf']['sum']['high'][0]-true_sum_highs_sf[j]) / x['sf']['sum']['high'][1] for x in a]) for j,a in enumerate(all_results)])
    pull_std_sf = np.asarray([np.std( [(x['sf']['sum']['high'][0]-true_sum_highs_sf[j]) / x['sf']['sum']['high'][1] for x in a]) for j,a in enumerate(all_results)])

    pull_means_of = np.asarray([np.mean( [(x['of']['sum']['high'][0]-true_sum_highs_of[j]) / x['of']['sum']['high'][1] for x in a]) for j,a in enumerate(all_results)])
    pull_std_of = np.asarray([np.std( [(x['of']['sum']['high'][0]-true_sum_highs_of[j]) / x['of']['sum']['high'][1] for x in a]) for j,a in enumerate(all_results)])

    means_of = np.asarray([np.median([x['of']['sum']['high'][0] for x in a]) for a in all_results])
    std_of = np.asarray([np.std([x['of']['sum']['high'][0] for x in a]) for a in all_results])

    pulls_sf = [
                [(a['sf']['sum']['high'][0]-true_sum_highs_sf[j])/a['sf']['sum']['high'][1] for a in all_results[j]]\
                for j in xrange(len(all_results))
               ]

    pulls_of = [
                [(a['of']['sum']['high'][0]-true_sum_highs_of[j])/a['of']['sum']['high'][1] for a in all_results[j]]\
                for j in xrange(len(all_results))
               ]

    # low = np.asarray([scipy.stats.scoreatpercentile([x['sf']['sum']['high'][0] for x in a], 16) for a in all_results])
    # high = np.asarray([scipy.stats.scoreatpercentile([x['sf']['sum']['high'][0] for x in a], 84) for a in all_results])


    plt.figure()
    plt.errorbar(true_sum_highs_sf, means_sf, yerr=std_sf, color="k", fmt='.')
    x = np.linspace(min(true_sum_highs_sf), max(true_sum_highs_sf), 500)
    plt.plot(x, x, '--', color='b')
    plt.legend(['Averaged Fitted Value (100 Toys)', 'Unit Slope'], loc="upper left")
    plt.xlabel("Generated Signal Events")
    plt.ylabel("Fitted Signal Events")
    plt.title("Same-Flavor")

    plt.figure()
    plt.errorbar(true_sum_highs_of, means_of, yerr=std_of, color="k", fmt='.')
    x = np.linspace(min(true_sum_highs_of), max(true_sum_highs_of), 500)
    plt.plot(x, x, '--', color='b')
    plt.legend(['Averaged Fitted Value (100 Toys)', 'Unit Slope'], loc="upper left")
    plt.xlabel("Generated Signal Events")
    plt.ylabel("Fitted Signal Events")
    plt.title("Opposite-Flavor")

    plt.figure()
    plt.hist(pulls_sf[0], bins=100, range=(-3,3))
    plt.xlabel("Signal Region Pull")
    plt.title("Same-Flavor")

    plt.figure()
    plt.hist(pulls_of[0], bins=100, range=(-3,3))
    plt.xlabel("Signal Region Pull")
    plt.title("Opposite-Flavor")
    plt.show()




def fit_toy(ws, n):

    high=200

    model = ws.obj("ModelConfig")

    ROOT.RooRandom.randomGenerator().SetSeed()

    ts = ROOT.RooStats.NumEventsTestStat(model.GetPdf())

    sampler = ROOT.RooStats.ToyMCSampler(ts, n)
    sampler.SetPdf(model.GetPdf())
    sampler.SetObservables(model.GetObservables())
    # sampler.SetGlobalObservables(model.GetGlobalObservables())
    sampler.SetParametersForTestStat(ROOT.RooArgSet())

    sampler.SetGenerateBinned(True)

    constr = model.GetNuisanceParameters()
    ROOT.RooStats.RemoveConstantParameters(constr)
    ws.obj("obs_x_sf").setRange("fit", 10, 120.)
    ws.obj("obs_x_of").setRange("fit", 10, 120.)

    results = []
    for _ in xrange(n):
        toy_data = sampler.GenerateToyData(model.GetSnapshot())
        res = model.GetPdf().fitTo(toy_data, ROOT.RooFit.Constrain(constr), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Range("fit"))

        r = ROOT.get_results(ws, res, 120., high)
        result = bkg_fit.results_to_dict_counting(r)
        result['of']['gen']['high'] = toy_data.sumEntries("(channelCat==channelCat::of) & (obs_x_of>120)")
        result['sf']['gen']['high'] = toy_data.sumEntries("(channelCat==channelCat::sf) & (obs_x_sf>120)")


        results.append(result)

    return results


if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)
    print args

    file_name = args['<filename>']



    # ncpu = int(args['--ncpu'])

    res = run_cut_based(file_name, 0)

