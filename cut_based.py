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

# rc = Client()
# dview = rc[:]
# with dview.sync_imports():
import my_lib
import ROOT
# dview.execute("ROOT.gROOT.ProcessLineSync('.L IntegralError.C+')")

ROOT.gROOT.ProcessLineSync('.L IntegralError.C+')

def run_cut_based(file_name, ncpu):

    rfile = ROOT.TFile(file_name)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    for x in xrange(11, 29):
        try:
            constr.find("gamma_z_syst_sf_bin_"+str(x)).setConstant(True)
        except AttributeError:
            pass
        try:
            constr.find("gamma_z_syst_of_bin_"+str(x)).setConstant(True)
        except AttributeError:
            pass
        try:
            constr.find("gamma_ww_syst_of_bin_"+str(x)).setConstant(True)
        except AttributeError:
            pass
        try:
            constr.find("gamma_ww_syst_sf_bin_"+str(x)).setConstant(True)
        except AttributeError:
            pass
        try:
            constr.find("gamma_wjets_syst_of_bin_"+str(x)).setConstant(True)
        except AttributeError:
            pass
        try:
            constr.find("gamma_wjets_syst_sf_bin_"+str(x)).setConstant(True)
        except AttributeError:
            pass
        try:
            constr.find("gamma_stat_of_bin_"+str(x)).setConstant(True)
        except AttributeError:
            pass
        try:
            constr.find("gamma_stat_sf_bin_"+str(x)).setConstant(True)
        except AttributeError:
            pass


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


    r = ROOT.get_results(ws, res)
    data_results = results_to_dict(r)
    print json.dumps(data_results, indent=4)

    n_events = model.GetPdf().getPdf("of").expectedEvents(ROOT.RooArgSet(obs_of))+\
                model.GetPdf().getPdf("sf").expectedEvents(ROOT.RooArgSet(obs_sf))

    n_toy_events = 10000
    true = data_results['sf']['sum']['high'][0]*n_toy_events/n_events

    temp_file = ROOT.TFile("temp.root", "recreate")
    ws.Write()
    temp_file.Close()
    # local test
    print fit_toy("temp.root", 1000000000)

    # lview = rc.load_balanced_view()

    means = []
    widths=[]

    # n = [10000, 50000, 70000, 100000, 200000, 300000, 400000, 500000, 700000, 1000000]
    n = [10000, 1000000]

    # all_results = defaultdict(list)

    # for n_toy_events in n:
    #     results = []
    #     true = data_results['sf']['sum']['high'][0]*n_toy_events/n_events
    #     # true = data_results['sf']['sum']['high'][0]


    #     for i in xrange(100):
    #         r = lview.apply_async(fit_toy, "temp.root", n_toy_events)
    #         results.append(r)

    #     lview.wait(results)

    #     sf_vals = []
    #     sf_errs = []
    #     for r in results:
    #             all_results[n_toy_events].append(r.result)
    #             sf_vals.append(r.result['sf']['sum']['high'][0])
    #             sf_errs.append(r.result['sf']['sum']['high'][1])

    #     dview.results.clear()
    #     # dview.clear(block=True)
    #     sf_vals = np.asarray(sf_vals)
    #     sf_errs = np.asarray(sf_errs)

    #     pulls = (sf_vals-true)/sf_vals.std()

    #     means.append(pulls.mean())
    #     widths.append(pulls.std())

    # means = np.asarray(means)
    # widths = np.asarray(widths)
    # plt.plot(n, means, color="k")
    # # plt.fill_between(n, means-widths, means+widths, color="b", alpha=0.5)
    # plt.show()


    # import IPython
    # IPython.embed()

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
    # ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


    # f = ROOT.TFile(ws_filename)
    # ws = f.Get("combined")
    # model = ws.obj("ModelConfig")

    # # ROOT.SetOwnership(ws, True)
    # ROOT.SetOwnership(model, True)

    # dummy = ROOT.RooStats.NumEventsTestStat(model.GetPdf()) # this doesn't do anything, we just need something to give ToyMCSampler
    # mc = ROOT.RooStats.ToyMCSampler(dummy, 1)
    # mc.SetPdf(model.GetPdf())
    # mc.SetObservables(model.GetObservables())
    # mc.SetGlobalObservables(model.GetGlobalObservables())
    # mc.SetParametersForTestStat(model.GetParametersOfInterest())
    # mc.SetNEventsPerToy(n)

    # constr = model.GetNuisanceParameters()
    # ROOT.RooStats.RemoveConstantParameters(constr)
    # ROOT.SetOwnership(constr, True)

    # toy_data = mc.GenerateToyData(model.GetSnapshot())
    # res = model.GetPdf().fitTo(toy_data, ROOT.RooFit.Constrain(constr), ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
    # ROOT.SetOwnership(toy_data, True)
    # ROOT.SetOwnership(res, True)
    
    # r = ROOT.get_results(ws, res)
    # result = my_lib.results_to_dict(r)
    # result['generated_sig']['sf'] = toy_data.sumEntries("(channelCat==channelCat::sf) & (obs_x_sf > 120.)")
    # result['generated_sig']['of'] = toy_data.sumEntries("(channelCat==channelCat::of) & (obs_x_of > 120.)")

    r = ROOT.fit_toy(ws_filename, n)
    result = my_lib.results_to_dict(r)
    return result

    # print "mc"
    # #mc.IsA().Destructor( mc )
    # print "dummy"
    # #dummy.IsA().Destructor( dummy )
    # print "constr"
    # constr.IsA().Destructor( constr )
    # print "toy_data"
    # toy_data.IsA().Destructor( toy_data )
    # print "res"
    # res.IsA().Destructor( res )
    # print "model"
    # model.IsA().Destructor( model )
    # print "ws"
    # ws.Delete()
    # print "done"

    # f.Close()

    # return result


if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)
    print args

    file_name = args['<filename>']



    # ncpu = int(args['--ncpu'])

    res = run_cut_based(file_name, 0)

