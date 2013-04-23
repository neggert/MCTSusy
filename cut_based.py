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

rc = Client()
dview = rc[:]
with dview.sync_imports():
    import my_lib
    import ROOT
dview.execute("ROOT.gROOT.ProcessLineSync('.L IntegralError.C+')")

ROOT.gROOT.ProcessLineSync('.L IntegralError.C+')

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

    obs_sf.setRange("fitRange_sf", 10., 120.)
    obs_of.setRange("fitRange_of", 10., 120.)


    # run the fit
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, ROOT.RooFit.Constrain(constr), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(0),
                               ROOT.RooFit.Range("fitRange"), ROOT.RooFit.SplitRange())

    params = model.GetParametersOfInterest()
    params.add(model.GetNuisanceParameters())
    model.SetSnapshot(params)

    temp_file = ROOT.TFile("temp.root", "recreate")
    ws.Write()
    temp_file.Close()

    r = ROOT.get_results(ws, res, 120)
    data_results = results_to_dict(r)
    print json.dumps(data_results, indent=4)

    # trial_cuts = np.arange(100, 300, 10)

    lview = rc.load_balanced_view()

    all_results = []

    # true_vals = [0 for _ in trial_cuts]

    # for j, trial_cut in enumerate(trial_cuts):
        # r = ROOT.get_results(ws, res, trial_cut)
        # true_vals[j] = results_to_dict(r)['sf']['sum']['high'][0]

    results = []
    
    for i in xrange(100):
        r = lview.apply_async(fit_toy, "temp.root", 0, 120.)
        results.append(r)

    lview.wait(results)

    for r in results:
            all_results.append(r.result)


    dview.results.clear()
        # dview.clear(block=True)

    means = np.mean([x['sf']['sum']['high'][0] for x in all_results])
    stds = np.std([x['sf']['sum']['high'][0] for x in all_results])

    low = scipy.stats.scoreatpercentile([x['sf']['sum']['high'][0] for x in all_results], 16)
    high = scipy.stats.scoreatpercentile([x['sf']['sum']['high'][0] for x in all_results], 84)

    # plt.plot(trial_cuts, means, color="k")
    # plt.fill_between(trial_cuts, low, high, color="b", alpha=0.5)
    # plt.plot(trial_cuts, true_vals, '--', color='k')
    # plt.show()


    import IPython
    IPython.embed()

def results_to_dict(r):
    results = defaultdict(lambda: defaultdict(dict))
    results['sf']['sum']['low'] = (r.sf.sum.low.val, r.sf.sum.low.error)
    results['sf']['sum']['high'] = (r.sf.sum.high.val, r.sf.sum.high.error)
    results['sf']['vv']['low'] = (r.sf.vv.low.val, r.sf.vv.low.error)
    results['sf']['vv']['high'] = (r.sf.vv.high.val, r.sf.vv.high.error)
    results['sf']['top']['low'] = (r.sf.top.low.val, r.sf.top.low.error)
    results['sf']['top']['high'] = (r.sf.top.high.val, r.sf.top.high.error)
    results['sf']['z']['low'] = (r.sf.z.low.val, r.sf.z.low.error)
    results['sf']['z']['high'] = (r.sf.z.high.val, r.sf.z.high.error)
    results['sf']['fake']['low'] = (r.sf.fake.low.val, r.sf.fake.low.error)
    results['sf']['fake']['high'] = (r.sf.fake.high.val, r.sf.fake.high.error)
    results['of']['sum']['low'] = (r.of.sum.low.val, r.of.sum.low.error)
    results['of']['sum']['high'] = (r.of.sum.high.val, r.of.sum.high.error)
    results['of']['vv']['low'] = (r.of.vv.low.val, r.of.vv.low.error)
    results['of']['vv']['high'] = (r.of.vv.high.val, r.of.vv.high.error)
    results['of']['top']['low'] = (r.of.top.low.val, r.of.top.low.error)
    results['of']['top']['high'] = (r.of.top.high.val, r.of.top.high.error)
    results['of']['z']['low'] = (r.of.z.low.val, r.of.z.low.error)
    results['of']['z']['high'] = (r.of.z.high.val, r.of.z.high.error)
    results['of']['fake']['low'] = (r.of.fake.low.val, r.of.fake.low.error)
    results['of']['fake']['high'] = (r.of.fake.high.val, r.of.fake.high.error)

    return results

def fit_toy(ws_filename, n, cut):

    r = ROOT.fit_toy(ws_filename, n, cut)
    result = my_lib.results_to_dict(r)
    return result


if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)
    print args

    file_name = args['<filename>']



    # ncpu = int(args['--ncpu'])

    res = run_cut_based(file_name, 0)

