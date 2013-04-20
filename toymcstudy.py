#!/usr/bin/env python

"""Set limits on the specified model

Usage:
    toymcstudy.py.py sig <model_file> <npoints> <tries_per_point>
    toymcstudy.py.py large <model_file> <tries_per_point>
    toymcstudy.py.py wjets <model_file> <npoints> <tries_per_point>
    toymcstudy.py.py vv <model_file> <tries_per_point>



Options:
    -h --help        Show this screen.

"""

from set_limits import *
from collections import defaultdict
import json
import numpy as np
import scipy.stats
import matplotlib.pylab as plt

def run_sig_pulls(model_file, npoints, tries_per_point):

    rfile = R.TFile(model_file)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    pars = model.GetNuisanceParameters()

    # run the initial fit
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.PrintLevel(0))

    poi_hat = model.GetParametersOfInterest().first().getVal()
    poi_hat_err = model.GetParametersOfInterest().first().getError()

    pois = np.linspace(0.001, 10, npoints)

    avg_list = []
    avg_err_up_list = []
    avg_err_down_list = []

    from IPython.parallel import Client

    rc = Client()
    dview = rc[:]

    dview.execute("import ROOT as R")
    lview = rc.load_balanced_view()

    for sig_strength in pois:

        results = []

        for i in xrange(50):
            r = lview.apply_async(get_sig_dist, ws, sig_strength, 10)
            results.append(r)

        lview.wait(results)

        sig_yields_fit = []
        for r in results:
            sig_yields_fit.extend(r.result)


        avg = np.median(sig_yields_fit)
        avg_err_up = scipy.stats.scoreatpercentile(sig_yields_fit, 84)
        avg_err_down = scipy.stats.scoreatpercentile(sig_yields_fit, 16)
        avg_list.append(avg)
        avg_err_up_list.append(avg_err_up)
        avg_err_down_list.append(avg_err_down)

    avg_list = np.asarray(avg_list)
    avg_err_up_list = np.asarray(avg_err_up_list)
    avg_err_down_list = np.asarray(avg_err_down_list)


    avg_list = np.asarray(avg_list)
    fig = plt.figure()
    plt.plot(pois, avg_list)
    plt.plot(pois, pois,'--')
    plt.fill_between(pois, avg_err_down_list, avg_err_up_list, alpha=0.3)
    plt.xlabel("Input Signal Strength")
    plt.ylabel("Fitted Signal Strength")
    plt.show()

    plt.savefig("plots/toymc.pdf")

    raw_input("...")


    # mc.fitParDataSet().Print()

def get_sig_dist(ws, poi, n):
    model = ws.obj("ModelConfig")

    params = model.GetParametersOfInterest()
    params.add(model.GetNuisanceParameters())
    model.GetParametersOfInterest().first().setVal(poi)
    model.SetSnapshot(params)

    dummy = R.RooStats.MinNLLTestStat(model.GetPdf()) # this doesn't do anything, we just need something to give ToyMCSampler
    mc = R.RooStats.ToyMCSampler(dummy, 1)
    mc.SetPdf(model.GetPdf())
    mc.SetObservables(model.GetObservables())
    mc.SetGlobalObservables(model.GetGlobalObservables())
    mc.SetParametersForTestStat(model.GetParametersOfInterest())

    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    sig_yields_fit = []
    for i in xrange(n):
        toy_data = mc.GenerateToyData(model.GetSnapshot())
        model.GetPdf().fitTo(toy_data, R.RooFit.Constrain(constr), R.RooFit.PrintLevel(-1))
        sig_yields_fit.append(model.GetParametersOfInterest().first().getVal())

    return sig_yields_fit

def run_large_sig(model_file, tries_per_point):

    rfile = R.TFile(model_file)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    pars = model.GetNuisanceParameters()

    # run the initial fit
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save())

    poi_hat = model.GetParametersOfInterest().first().getVal()
    poi_hat_err = model.GetParametersOfInterest().first().getError()


    params = model.GetParametersOfInterest()
    params.add(model.GetNuisanceParameters())
    model.GetParametersOfInterest().first().setVal(10.)
    model.SetSnapshot(params)

    sig_strength = 10.
    # now generate toys
    dummy = R.RooStats.MinNLLTestStat(model.GetPdf()) # this doesn't do anything, we just need something to give ToyMCSampler
    mc = R.RooStats.ToyMCSampler(dummy, 1)
    mc.SetPdf(model.GetPdf())
    mc.SetObservables(model.GetObservables())
    mc.SetGlobalObservables(model.GetGlobalObservables())
    mc.SetParametersForTestStat(model.GetParametersOfInterest())

    sig_yields_fit = []
    sig_yields_err_fit = []
    for i in xrange(tries_per_point):
        model.LoadSnapshot()
        toy_data = mc.GenerateToyData(params)
        model.GetPdf().fitTo(toy_data, R.RooFit.Constrain(constr), R.RooFit.PrintLevel(-1))
        sig_yields_fit.append(model.GetParametersOfInterest().first().getVal())
        sig_yields_err_fit.append(model.GetParametersOfInterest().first().getError())

    sig_yields_fit = np.asarray(sig_yields_fit)
    sig_yields_err_fit = np.asarray(sig_yields_err_fit)

    pulls = (sig_yields_fit-sig_strength)/sig_yields_err_fit

    print sig_yields_fit
    print sig_yields_err_fit


    avg = np.mean(sig_yields_fit)
    avg_err = np.std(sig_yields_fit)

    # plot pull
    plt.figure()
    plt.hist(pulls, range=(-3, 3), bins=30, histtype="stepfilled")
    plt.xlabel("Pull")

    # plot NLL for last toy
    nll = model.GetPdf().createNLL(toy_data, R.RooFit.Constrain(constr))
    poi = model.GetParametersOfInterest().first()
    f = poi.frame(sig_strength/10, sig_strength*10)
    nll.plotOn(f)
    f.Draw()

    import IPython
    IPython.embed()

def run_wjets_pulls(model_file, npoints, tries_per_point):

    rfile = R.TFile(model_file)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")
    model.GetParametersOfInterest().first().setVal(0.)
    model.GetParametersOfInterest().first().setConstant(True)

    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    pars = model.GetNuisanceParameters()

    err_pars = R.RooArgSet(pars.find("n_of_top"), pars.find("n_of_vv"), pars.find("n_of_z"), pars.find("n_of_wjets"),
                       pars.find("n_sf_top"), pars.find("n_sf_vv"), pars.find("n_sf_z"), pars.find("n_sf_wjets"))

    # run the initial fit
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save())

    wjets_hat = model.GetNuisanceParameters().find("n_of_wjets").getVal()
    wjets_hat_err = model.GetNuisanceParameters().find("n_of_wjets").getError()

    wjets_vals = np.linspace(max([0., wjets_hat-2*wjets_hat_err]), wjets_hat+5*wjets_hat_err, npoints)

    avg_list = []
    avg_err_up_list = []
    avg_err_down_list = []

    for wjets_val in wjets_vals:

        params = model.GetParametersOfInterest()
        params.add(model.GetNuisanceParameters())
        model.GetNuisanceParameters().find('n_of_wjets').setVal(wjets_val)
        model.SetSnapshot(params)


        # now generate toys
        dummy = R.RooStats.MinNLLTestStat(model.GetPdf()) # this doesn't do anything, we just need something to give ToyMCSampler
        mc = R.RooStats.ToyMCSampler(dummy, 1)
        mc.SetPdf(model.GetPdf())
        mc.SetObservables(model.GetObservables())
        mc.SetGlobalObservables(model.GetGlobalObservables())
        mc.SetParametersForTestStat(model.GetParametersOfInterest())

        yields_fit = []
        yields_err_fit = []
        for i in xrange(tries_per_point):
            model.LoadSnapshot()
            toy_data = mc.GenerateToyData(params)
            model.GetPdf().fitTo(toy_data, R.RooFit.Constrain(constr), R.RooFit.PrintLevel(-1))
            yields_fit.append(model.GetNuisanceParameters().find('n_of_wjets').getVal())
            yields_err_fit.append(model.GetNuisanceParameters().find('n_of_wjets').getError())

        yields_fit = np.asarray(yields_fit)
        yields_err_fit = np.asarray(yields_err_fit)

        print yields_fit
        print yields_err_fit

        avg = np.median(yields_fit)
        avg_err_up = scipy.stats.scoreatpercentile(yields_fit, 84)
        avg_err_down = scipy.stats.scoreatpercentile(yields_fit, 16)

        avg_list.append(avg)
        avg_err_up_list.append(avg_err_up)
        avg_err_down_list.append(avg_err_down)

    avg_list = np.asarray(avg_list)
    avg_err_up_list = np.asarray(avg_err_up_list)
    avg_err_down_list = np.asarray(avg_err_down_list)

    fig = plt.figure()
    plt.plot(wjets_vals, avg_list)
    plt.plot(wjets_vals, wjets_vals,'--')
    plt.fill_between(wjets_vals, avg_err_down_list, avg_err_up_list, alpha=0.3)
    plt.xlabel("Input WJets yield")
    plt.ylabel("Fitted Wjets yield")
    plt.show()

    plt.savefig("plots/toymc_wjets.pdf")

    raw_input("...")

def run_sfvv_pulls(model_file, tries_per_point):

    rfile = R.TFile(model_file)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")
    model.GetParametersOfInterest().first().setVal(0.)
    model.GetParametersOfInterest().first().setConstant(True)

    constr = model.GetNuisanceParameters()
    n_sf_vv = constr.find("n_vv_sf")

    R.RooStats.RemoveConstantParameters(constr)



    # run the initial fit
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.PrintLevel(0))


    vv_list = []
    n_high_list = []
    toys = []


    params = model.GetParametersOfInterest()
    params.add(model.GetNuisanceParameters())
    model.SetSnapshot(params)
    
    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    obs_sf = ws.obj("obs_x_sf")
    obs_of = ws.obj("obs_x_of")



    # now generate toys
    dummy = R.RooStats.MinNLLTestStat(model.GetPdf()) # this doesn't do anything, we just need something to give ToyMCSampler
    mc = R.RooStats.ToyMCSampler(dummy, 1)
    mc.SetPdf(model.GetPdf())
    mc.SetObservables(model.GetObservables())
    mc.SetGlobalObservables(model.GetGlobalObservables())
    mc.SetParametersForTestStat(model.GetParametersOfInterest())

    for i in xrange(tries_per_point):
        model.LoadSnapshot()
        toy_data = mc.GenerateToyData(params)
        toys.append(toy_data)
        model.GetPdf().fitTo(toy_data, R.RooFit.Constrain(constr), R.RooFit.PrintLevel(-1))
        vv_list.append(n_sf_vv.getVal())
        n_high = toy_data.createHistogram(obs_sf, obs_of).ProjectionX("t_of"+str(i), 28,29).Integral(13, 29)
        n_high_list.append(n_high)

    fig = plt.figure()
    plt.hist(vv_list, histtype="stepfilled", range=(0, 3000), bins=30)
    plt.xlabel("SF Diboson yield")
    plt.show()

    plt.savefig("plots/toymc_sfvv.pdf")


    # munge this to get it into the correct form for boxplot
    fitted_values_vs_nhigh = defaultdict(list)
    for i in xrange(tries_per_point):
        fitted_values_vs_nhigh[n_high_list[i]].append(vv_list[i])
    fig = plt.figure()
    plt.boxplot([fitted_values_vs_nhigh[key] for key in sorted(fitted_values_vs_nhigh.iterkeys())])
    plt.ylabel("SF Diboson yield")
    plt.xlabel("Events w/ MCTPerp > 130 GeV")
    plt.show()

    plt.savefig("plots/toymc_sfvv_vs_nhigh.pdf")


if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)

    try:
        npoints = int(args['<npoints>'])
    except TypeError:
        pass
    tries_per_point = int(args['<tries_per_point>'])

    model_file = args['<model_file>']

    if args['sig']:
        run_sig_pulls(model_file, npoints, tries_per_point)
    elif args['large']:
        run_large_sig(model_file, tries_per_point)
    elif args['wjets']:
        run_wjets_pulls(model_file, npoints, tries_per_point)
    elif args['vv']:
        run_sfvv_pulls(model_file, tries_per_point)
