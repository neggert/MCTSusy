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
    params = model.GetParametersOfInterest()
    params.add(model.GetNuisanceParameters())
    model.SetSnapshot(params)

    # run the fit
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, ROOT.RooFit.Constrain(constr), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(0),
                               ROOT.RooFit.Range("fitRange"), ROOT.RooFit.SplitRange())



    r = ROOT.get_results(ws, res)
    data_results = results_to_dict(r)
    print json.dumps(data_results, indent=4)

    # n_sf_vv_nominal = ws.obj("n_vv_sf").getVal()
    # n_sf_top_nominal = ws.obj("n_top_sf").getVal()

    lview = rc.load_balanced_view()

    n_tries = 20
    all_results = [[] for _ in xrange(n_tries)]

    true_sum_highs = []

    for j in xrange(n_tries):

        model.LoadSnapshot()
        rv = generate_random_norms(5393,5314)
        print rv
        ws.obj("n_top_sf").setVal(rv['n_top_sf'])
        ws.obj("n_sf_z").setVal(rv['n_z_sf'])
        ws.obj("n_sf_wjets").setVal(rv['n_fake_sf'])
        ws.obj("n_of_z").setVal(rv['n_z_of'])
        ws.obj("n_of_wjets").setVal(rv['n_fake_of'])
        ws.obj("alpha_t_vv_ratio_sf").setVal(0.)
        ws.obj("alpha_top_ratio").setVal(0.)
        ws.obj("alpha_vv_ratio").setVal(0.)
        model.SetSnapshot(params)

        r = ROOT.get_results(ws, res)
        data_results = results_to_dict(r)
        true_sum_highs.append(data_results['sf']['sum']['high'][0])

        model.SetSnapshot(params)
        temp_file = ROOT.TFile("temp.root", "recreate")
        ws.Write()
        temp_file.Close()

        results = []
        
        for i in xrange(100):
            r = lview.apply_async(fit_toy, "temp.root", 0)
            results.append(r)

        lview.wait(results)

        for r in results:
                all_results[j].append(r.result)


        dview.results.clear()
        # dview.clear(block=True)

    true_sum_highs = np.asarray(true_sum_highs)

    means = np.asarray([np.median([x['sf']['sum']['high'][0] for x in a]) for a in all_results])
    std = np.asarray([np.median([x['sf']['sum']['high'][0] for x in a]) for a in all_results])/np.sqrt(100)

    low = np.asarray([scipy.stats.scoreatpercentile([x['sf']['sum']['high'][0] for x in a], 16) for a in all_results])
    high = np.asarray([scipy.stats.scoreatpercentile([x['sf']['sum']['high'][0] for x in a], 84) for a in all_results])

    sort_order = np.argsort(true_sum_highs)
    true_sum_highs = true_sum_highs[sort_order]
    means = means[sort_order]
    low = low[sort_order]
    high = high[sort_order]

    plt.errorbar(true_sum_highs, means, yerr=std, color="k")
    # plt.fill_between(true_sum_highs, low, high, color="b", alpha=0.5)
    x = np.linspace(min(true_sum_highs), max(true_sum_highs), 500)
    plt.plot(x, x, '--', color='k')
    plt.show()


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

def fit_toy(ws_filename, n):

    r = ROOT.fit_toy(ws_filename, n)
    result = my_lib.results_to_dict(r)
    return result

def generate_random_norms(n_sf, n_of):
    n_top_sf = np.random.random()*n_sf
    n_top_of = n_top_sf*1.35275
    n_vv_sf = n_top_sf*0.502
    n_vv_of = n_vv_sf*1.06567

    n_sf_left = n_sf - n_top_sf - n_vv_sf
    n_of_left = n_of - n_top_of - n_vv_of

    if n_sf_left < 0 or n_of_left < 0:
        return generate_random_norms(n_sf, n_of)

    n_z_sf = np.random.random()*n_sf_left
    n_fake_sf = n_sf_left-n_z_sf

    n_z_of = np.random.random()*n_of_left
    n_fake_of = n_of_left - n_z_of

    return {"n_top_sf":n_top_sf,
            "n_top_of": n_top_of,
            "n_vv_sf":n_vv_sf,
            "n_vv_of":n_vv_of,
            "n_z_sf":n_z_sf,
            "n_z_of":n_z_of,
            "n_fake_of":n_fake_of,
            "n_fake_sf":n_fake_sf
    }

if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)
    print args

    file_name = args['<filename>']



    # ncpu = int(args['--ncpu'])

    res = run_cut_based(file_name, 0)

