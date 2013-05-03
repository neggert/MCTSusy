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

import my_lib
import ROOT

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
    # constrPdf = ws.obj("constrPdf")

    constr = model.GetNuisanceParameters()

    top_constraint = ws.obj("top_ratio_constraint")
    vv_constraint = ws.obj("vv_ratio_constraint")
    top_vv_constraint_sf = ws.obj("top_vv_ratio_sf_constraint")
    top_vv_constraint_of = ws.obj("top_vv_ratio_of_constraint")

    ex_const = ROOT.RooArgSet(top_constraint, vv_constraint, top_vv_constraint_sf, top_vv_constraint_of)



    ROOT.RooStats.RemoveConstantParameters(constr)


    model.GetParametersOfInterest().first().setVal(0.)
    model.GetParametersOfInterest().first().setConstant()

    obs_sf = ws.obj("obs_x_sf")
    obs_of = ws.obj("obs_x_of")

    obs_sf.setRange("fitRange", 10., 120.)
    obs_of.setRange("fitRange", 10., 120.)

    params = model.GetParametersOfInterest()
    params.add(model.GetNuisanceParameters())


    # run the fit
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(0), ROOT.RooFit.ExternalConstraints(ex_const),
                               ROOT.RooFit.Range("fitRange"), ROOT.RooFit.Constrain(constr))
    model.SetSnapshot(params)



    r = ROOT.get_results(ws, res)
    data_results = results_to_dict(r)
    # print json.dumps(data_results, indent=4)

    # test
    # fit_toy(ws, 0)

    # n_sf_vv_nominal = ws.obj("n_vv_sf").getVal()
    # n_sf_top_nominal = ws.obj("n_top_sf").getVal()

    lview = rc.load_balanced_view()

    n_tries = 1
    all_results = [[] for _ in xrange(n_tries)]

    true_sum_highs_sf = []
    true_sum_highs_of = []





    for j in xrange(n_tries):

        model.LoadSnapshot()
        # rv = generate_random_norms(5393,5314)
        # rv = {'n_top_sf':2271,
        #       'n_top_of':3071,
        #       'n_vv_sf':1249,
        #       'n_vv_of':1466,
        #      'n_z_sf':920,
        #      'n_fake_sf': 607,
        #      'n_z_of': 58,
        #      'n_fake_of': 0.1
        #      }
        # print rv
        # ws.obj("n_sf_top").setVal(rv['n_top_sf'])
        # ws.obj("n_of_top").setVal(rv['n_top_of'])
        # ws.obj("n_sf_vv").setVal(rv['n_vv_sf'])
        # ws.obj("n_of_vv").setVal(rv['n_vv_of'])
        # ws.obj("n_sf_z").setVal(rv['n_z_sf'])
        # ws.obj("n_sf_wjets").setVal(rv['n_fake_sf'])
        # ws.obj("n_of_z").setVal(rv['n_z_of'])
        # ws.obj("n_of_wjets").setVal(rv['n_fake_of'])
        model.SetSnapshot(params)

        r = ROOT.get_results(ws, res)
        data_results = results_to_dict(r)
        true_sum_highs_sf.append(data_results['sf']['sum']['high'][0])
        true_sum_highs_of.append(data_results['of']['sum']['high'][0])

        results = []
        
        for i in xrange(1000):
            r = lview.apply_async(fit_toy, ws, 0)
            results.append(r)

        lview.wait(results)

        for r in results:
                all_results[j].append(r.result)


        dview.results.clear()
        # dview.clear(block=True)

    true_sum_highs_sf = np.asarray(true_sum_highs_sf)
    true_sum_highs_of = np.asarray(true_sum_highs_of)



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

def fit_toy(ws, n):

    # globals = ROOT.RooArgSet()
    # for i in xrange(29):
    #     globals.add(ws.obj("nom_gamma_z_syst_sf_bin_{}".format(i)))
    #     globals.add(ws.obj("nom_gamma_wjets_syst_sf_bin_{}".format(i)))
    #     globals.add(ws.obj("nom_gamma_stat_sf_bin_{}".format(i)))
    # for i in xrange(10):
    #     globals.add(ws.obj("nom_gamma_ww_syst_sf_bin_{}".format(i)))
    # globals.add(ws.obj("nom_alpha_HWW_norm"))
    # globals.add(ws.obj("nom_alpha_VVV_norm"))
    # globals.add(ws.obj("nom_alpha_WW_norm"))
    # globals.add(ws.obj("nom_alpha_WZ_norm"))
    # globals.add(ws.obj("nom_alpha_ZZ_norm"))
    # globals.add(ws.obj("nom_alpha_t_vv_ratio_sf"))
    # globals.add(ws.obj("nom_alpha_top_ratio"))
    # globals.add(ws.obj("nom_alpha_vv_ratio"))

    globals = ws.obj("ModelConfig").GetGlobalObservables()

    r = ROOT.fit_toy(ws, n, globals)
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

