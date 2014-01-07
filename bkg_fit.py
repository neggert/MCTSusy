#!/usr/bin/env python

#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    bkg_fit.py <filename> <outfile> [-pm] [--ncpu=<c>] [--channels=<c1,c2>] [--cut=<cut>]

Options:
    -h --help        Show this screen.
    --ncpu=<c>       Number of CPUs to use [default: 1]
    --channels=<c1,c2> Channels to use [default: of,sf]
    --cut=<cut>      MCTPerp value specifying signal region in counting analysis [default: -1]
    -p               Get p value (slow)
    -m               Run MINOS (slow)

"""

from set_limits import *
from collections import defaultdict
import json
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import IndexLocator, FixedFormatter
import numpy as np
from prettytable import PrettyTable
import money_take2
import ROOT
import sys
# plt.switch_backend("pdf")

bkgs = ['top', 'vv', 'wjets', 'z']

ROOT.gROOT.ProcessLineSync('.L IntegralErrorShape.C+')

def get_step_fill_between(x, y1, y2):

    fill_x = np.zeros(2*len(x))
    fill_x[::2] = x
    fill_x[1:-1:2] = x[1:]
    fill_x[-1] = 300
    fill_y1 = np.zeros(fill_x.shape)
    fill_y1[::2] = y1
    fill_y1[1::2] = y1
    fill_y2 = np.zeros(fill_x.shape)
    fill_y2[::2] = y2
    fill_y2[1::2] = y2



    return np.append(fill_x, fill_x[::-1]), np.append(fill_y1, fill_y2[::-1])

def results_to_dict(r):
    results = defaultdict(lambda: defaultdict(dict))
    results['sf']['sum'] = (r.sf.sum.low.val, r.sf.sum.low.error)
    results['sf']['vv'] = (r.sf.vv.low.val, r.sf.vv.low.error)
    results['sf']['top'] = (r.sf.top.low.val, r.sf.top.low.error)
    results['sf']['z'] = (r.sf.z.low.val, r.sf.z.low.error)
    results['sf']['fake'] = (r.sf.fake.low.val, r.sf.fake.low.error)
    results['of']['sum'] = (r.of.sum.low.val, r.of.sum.low.error)
    results['of']['vv'] = (r.of.vv.low.val, r.of.vv.low.error)
    results['of']['top'] = (r.of.top.low.val, r.of.top.low.error)
    results['of']['z'] = (r.of.z.low.val, r.of.z.low.error)
    results['of']['fake'] = (r.of.fake.low.val, r.of.fake.low.error)

    return results

def results_to_dict_counting(r):
    results = defaultdict(lambda: defaultdict(dict))
    results['sf']['sum']['low'] = (r.sf.sum.low.val, r.sf.sum.low.error)
    results['sf']['vv']['low'] = (r.sf.vv.low.val, r.sf.vv.low.error)
    results['sf']['top']['low'] = (r.sf.top.low.val, r.sf.top.low.error)
    results['sf']['z']['low'] = (r.sf.z.low.val, r.sf.z.low.error)
    results['sf']['fake']['low'] = (r.sf.fake.low.val, r.sf.fake.low.error)
    results['of']['sum']['low'] = (r.of.sum.low.val, r.of.sum.low.error)
    results['of']['vv']['low'] = (r.of.vv.low.val, r.of.vv.low.error)
    results['of']['top']['low'] = (r.of.top.low.val, r.of.top.low.error)
    results['of']['z']['low'] = (r.of.z.low.val, r.of.z.low.error)
    results['of']['fake']['low'] = (r.of.fake.low.val, r.of.fake.low.error)

    results['sf']['sum']['high'] = (r.sf.sum.high.val, r.sf.sum.high.error)
    results['sf']['vv']['high'] = (r.sf.vv.high.val, r.sf.vv.high.error)
    results['sf']['top']['high'] = (r.sf.top.high.val, r.sf.top.high.error)
    results['sf']['z']['high'] = (r.sf.z.high.val, r.sf.z.high.error)
    results['sf']['fake']['high'] = (r.sf.fake.high.val, r.sf.fake.high.error)
    results['of']['sum']['high'] = (r.of.sum.high.val, r.of.sum.high.error)
    results['of']['vv']['high'] = (r.of.vv.high.val, r.of.vv.high.error)
    results['of']['top']['high'] = (r.of.top.high.val, r.of.top.high.error)
    results['of']['z']['high'] = (r.of.z.high.val, r.of.z.high.error)
    results['of']['fake']['high'] = (r.of.fake.high.val, r.of.fake.high.error)

    return results

def run_bonly_fit(file_name, out_file, ncpu, get_p, data_prefix="data", data_file_name="data.root", do_minos=False, cut=None):

    high = 200.
    rfile = R.TFile(file_name)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    # model.GetParametersOfInterest().first().setVal(0.)
    # model.GetParametersOfInterest().first().setConstant()

    pars = model.GetNuisanceParameters()

    # err_pars = R.RooArgSet(pars.find("n_of_top"), pars.find("n_of_vv"), pars.find("n_of_z"), pars.find("n_of_wjets"),
    #                    pars.find("n_of_top"), pars.find("n_of_vv"), pars.find("n_of_z"), pars.find("n_of_wjets"))

    # run the fit
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR)
    counting = False
    if cut:
        ws.obj("obs_x_sf").setRange("fit", 10, cut)
        ws.obj("obs_x_of").setRange("fit", 10, cut)
        counting = True
    else:
        cut = high
    if do_minos:
        res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.PrintLevel(0), R.RooFit.Minos(), R.RooFit.Hesse(), R.RooFit.Range("fit"))
    else: 
        res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.PrintLevel(0), R.RooFit.Range("fit"))

    r = ROOT.get_results(ws, res, cut, high)
    # import pdb
    # pdb.set_trace()
    data_results = results_to_dict_counting(r) if counting else results_to_dict(r)


    fitPars = res.floatParsFinal()

    nPar = fitPars.getSize()

    t = PrettyTable()
    if do_minos:
        t.field_names = ['#', '', "Value", "Parabolic Error", "Minos Down", "Minos Up"]
    else:
        t.field_names = ['#', '', "Value", "Parabolic Error"]
    t.vertical_char = "&"
    for i in xrange(nPar):
        name = fitPars.at(i).GetName()
        val =  fitPars.at(i).getVal()
        err =  fitPars.at(i).getError()
        if do_minos:
            minos_up = fitPars.at(i).getErrorLo()
            minos_down = fitPars.at(i).getErrorHi()
            t.add_row((i, name, "{0:.3g}".format(val), "{0:.3g}".format(err), "{0:.3g}".format(minos_down), "{0:.3g}".format(minos_up)))

        else:
            t.add_row((i, name, "{0:.3g}".format(val), "{0:.3g}".format(err)))

    print t

    #find the actual parameter vlaues
    # par = R.RooArgList(ws.obj("n_top_sf"), ws.obj("alpha_top_ratio"))
    # n_top_sf_real = R.RooFormulaVar("n_top_sf_real", "n_top_sf_real", "n_top_sf*(1-0.1*alpha_top_ratio)", par)
    # scale = ws.obj("n_top_of_scale").getVal()
    # n_top_of = R.RooFormulaVar("n_top_of", "n_top_of", "n_top_sf*{0}*(1+0.1*alpha_top_ratio)".format(scale), par)

    # par = R.RooArgList(ws.obj("n_vv_sf"), ws.obj("alpha_vv_ratio"))
    # vv_ratio = ws.obj("alpha_vv_ratio")
    # n_vv_sf = ws.obj("n_vv_sf")
    # n_vv_sf_real = R.RooFormulaVar("n_vv_sf_real", "n_vv_sf_real", "n_vv_sf*(1-0.1*alpha_vv_ratio)", par)
    # n_vv_sf_err = n_vv_sf_real.getVal()*np.sqrt((n_vv_sf.getError()/n_vv_sf.getVal())**2+(0.1*vv_ratio.getError()/(1-0.1*vv_ratio.getVal()))**2)

    # scale = ws.obj("n_vv_of_scale").getVal()
    # n_vv_of = R.RooFormulaVar("n_vv_of", "n_vv_of", "n_vv_sf*{0}*(1+0.1*alpha_vv_ratio)".format(scale), par)

    # n_vv_of_err = n_vv_of.getVal()*np.sqrt((n_vv_sf.getError()/n_vv_sf.getVal())**2+(0.1*vv_ratio.getError()/(1+0.1*vv_ratio.getVal()))**2)

    f = open(out_file, 'w')

    json.dump(data_results, f, indent=3)

    f.close()

    # # plot the relevant portion of the correlation matrix
    # fullcor = res.correlationMatrix()
    # cor = fullcor.GetSub(129, 136, 129, 136)
    # # import pdb; pdb.set_trace()
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.patch.set_facecolor('gray')
    # ax.set_aspect('equal', 'box')
    # labels = FixedFormatter(['',
    #                         'OF Top',
    #                          'OF Diboson',
    #                          'OF Z',
    #                          'OF WJets',
    #                          'SF Top',
    #                          'SF Diboson',
    #                          'SF Z',
    #                          'SF WJets'])
    # ax.xaxis.set_major_formatter(labels)
    # ax.xaxis.tick_top()
    # ax.yaxis.set_major_formatter(labels)
    # for t in ax.get_xticklabels():
    #     t.set_rotation(40)
    #     t.set_ha('left')


    # # make Hinton-style correlation plot
    # for i in xrange(cor.GetNrows()):
    #     for j in xrange(cor.GetNcols()):
    #         # if i<=j: continue
    #         c = cor[i][j]
    #         if abs(c) < 0.01: continue
    #         if c > 0: color='white'
    #         else: color='black'
    #         size = np.sqrt(np.abs(c))
    #         rect = Rectangle([i-size/2, j-size/2], size, size, facecolor=color, edgecolor='black', lw=0.1)
    #         ax.add_patch(rect)
    # ax.autoscale_view()
    # plt.gca().invert_yaxis()
    # plt.tight_layout()

    # plt.savefig("plots/correlation.pdf")

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.patch.set_facecolor('gray')
    # ax.set_aspect('equal', 'box')

    # cor.Print()


    # # make Hinton-style correlation plot
    # for i in xrange(fullcor.GetNrows()):
    #     for j in xrange(fullcor.GetNcols()):
    #         # if i<=j: continue
    #         c = fullcor[i][j]
    #         if abs(c) < 0.01: continue
    #         if c > 0: color='white'
    #         else: color='black'
    #         size = np.sqrt(np.abs(c))
    #         rect = Rectangle([i-size/2, j-size/2], size, size, facecolor=color, edgecolor='black', lw=0.1)
    #         ax.add_patch(rect)
    # ax.set_xlim(0, fullcor.GetNrows())
    # ax.set_ylim(0, fullcor.GetNrows())
    # plt.gca().invert_yaxis()
    # plt.tight_layout()

    # plt.savefig("plots/correlation_full.pdf")
    # # raw_input("...")

    params = R.RooArgSet()
    params.add(model.GetNuisanceParameters())
    params.add(model.GetParametersOfInterest())
    model.SetSnapshot(params)

    money_take2.build_background_shape(ws, 'sf', money_take2.sf_backgrounds, log=True)
    money_take2.build_background_shape(ws, 'of', money_take2.of_backgrounds, log=True)
    money_take2.build_background_shape(ws, 'sf', money_take2.sf_backgrounds, log=False)
    money_take2.build_background_shape(ws, 'of', money_take2.of_backgrounds, log=False)


    # plot_fitted_sf(ws)
    # plot_fitted_of(ws)


    # # plot the fitted templates

    # # get the test statistic on data    
    R.gROOT.ProcessLineSync(".L KS/AndersonDarlingTestStat.cc+")
    AD = R.RooStats.AndersonDarlingTestStat(model.GetPdf())
    ts = AD.Evaluate(data)

    # # import IPython
    # # IPython.embed()


    # calculate a p-value
    if get_p:

        results = []
        from IPython.parallel import Client

        rc = Client()
        dview = rc[:]
        # with dview.sync_imports(): 
        #     import ROOT as R
        dview.execute("import ROOT")
        dview.execute("ROOT.gROOT.ProcessLineSync('.L KS/AndersonDarlingTestStat.cc+')")
        lview = rc.load_balanced_view()

        get_p_value_dist(ws, 2)

        for i in xrange(50):
            r = lview.apply_async(get_p_value_dist, ws, 100)
            results.append(r)

        lview.wait(results)

        sampDist = R.RooStats.SamplingDistribution()
        for r in results:
            sampDist.Add(r.result)


        f = ROOT.TFile("BkgADDist_test.root", "RECREATE")
        sampDist.Write("sampDist")
        f.Close()

        p = 1-sampDist.CDF(ts)

        print "P value:", p

    print "Test statistic on data: {:.7f}".format(ts)

    return data_results

def get_p_value_dist(ws, n):
    model = ws.obj("ModelConfig")

    ROOT.RooRandom.randomGenerator().SetSeed()

    AD = ROOT.RooStats.AndersonDarlingTestStat(model.GetPdf())

    sampler = ROOT.RooStats.ToyMCSampler(AD, n)
    sampler.SetPdf(model.GetPdf())
    sampler.SetObservables(model.GetObservables())
    sampler.SetGlobalObservables(model.GetGlobalObservables())
    sampler.SetParametersForTestStat(ROOT.RooArgSet())

    sampDist = sampler.GetSamplingDistribution(model.GetSnapshot())

    return sampDist

def plot_fitted_sf(ws):
    obs = ws.obj("obs_x_sf")
    x = np.arange(10, 300, 10)


    # statistical uncertainties
    stat_uncertainty_means = np.zeros(x.shape)
    for i in xrange(len(x)):
        try:
            stat_uncertainty_means[i] = ws.obj("gamma_stat_sf_bin_{0}_tau".format(i)).getVal()
        except AttributeError:
            pass

    stat_uncertainty_rel = 1./np.sqrt(stat_uncertainty_means)

    # top
    par = R.RooArgList(ws.obj("n_top_sf"), ws.obj("alpha_top_ratio"))
    n_top_sf_real = R.RooFormulaVar("n_top_sf_real", "n_top_sf_real", "n_top_sf*(1-0.1*alpha_top_ratio)", par)

    n_sf_top = n_top_sf_real.getVal()
    top_fitted = ws.obj("top_sf_sf_overallSyst_x_StatUncert")
    top_nominal = ws.obj("top_sf_sf_nominal")

    # import IPython
    # IPython.embed()

    top_fitted_points = np.zeros(x.shape)
    top_nom_points = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        top_fitted_points[i] = top_fitted.getVal()
        top_nom_points[i] = top_nominal.getVal()

    top_nom_points *= n_sf_top
    top_stat_high = top_nom_points*(1+stat_uncertainty_rel)
    top_stat_low = top_nom_points*(1-stat_uncertainty_rel)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill(*get_step_fill_between(x, top_stat_low, top_stat_high), color='r', alpha=0.5)
    ax.step(x, top_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, top_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Top SF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Statistical Systematic"])

    plt.savefig("plots/template_top_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill(*get_step_fill_between(x, top_stat_low, top_stat_high), color='r', alpha=0.5)
    ax.step(x, top_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, top_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_ylim(0.001, 5000)
    ax.set_title("Top SF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Statistical Systematic"])

    plt.savefig("plots/template_top_sf_log.pdf")

    # VV

    # fitted shape
    # includes both fitted statistical and systematics
    # might need to multiply by bin width
    par = R.RooArgList(ws.obj("n_vv_sf"), ws.obj("alpha_vv_ratio"))
    n_vv_sf_real = R.RooFormulaVar("n_vv_sf_real", "n_vv_sf_real", "n_vv_sf*(1-0.1*alpha_vv_ratio)", par)
    n_sf_vv = n_vv_sf_real.getVal()
    vv_fitted = ws.obj("vv_sf_sf_overallSyst_x_StatUncert_x_sf_ww_syst_sf_ShapeSys")
    vv_syst_nom = ws.obj("vv_sf_sf_Hist_alphanominal")

    n_vv_syst = 4
    vv_systs_low = []
    vv_systs_high = []

    for i in xrange(n_vv_syst):
        vv_systs_low.append(ws.obj("vv_sf_sf_Hist_alpha_{0}low".format(i)))
        vv_systs_high.append(ws.obj("vv_sf_sf_Hist_alpha_{0}high".format(i)))


    vv_fitted_points = np.zeros(x.shape)
    vv_nom_points = np.zeros(x.shape)

    vv_shape_syst_points = np.zeros(x.shape)

    vv_syst_low_points = np.zeros((n_vv_syst, len(x)))
    vv_syst_high_points = np.zeros((n_vv_syst, len(x)))

    for i in xrange(len(x)):
        obs.setVal(x[i])
        vv_fitted_points[i] = vv_fitted.getVal()
        vv_nom_points[i] = vv_syst_nom.getVal()
        try:
            vv_shape_syst_points[i] = ws.obj("gamma_ww_syst_sf_bin_{0}_sigma".format(i)).getVal()
        except AttributeError:
            pass
        for j in xrange(n_vv_syst):
            vv_syst_low_points[j,i] = vv_systs_low[j].getVal()
            vv_syst_high_points[j,i] = vv_systs_high[j].getVal()


    # now get the deviations from the nominal
    vv_syst_dev_low = vv_syst_low_points - vv_nom_points
    vv_syst_dev_high = vv_syst_high_points - vv_nom_points

    vv_total_syst_dev_low = np.sqrt(np.sum(vv_syst_dev_low**2, axis=0))
    vv_total_syst_dev_high = np.sqrt(np.sum(vv_syst_dev_high**2, axis=0))


    vv_syst_low_points =  (vv_nom_points - vv_total_syst_dev_low)*n_sf_vv
    vv_syst_high_points = (vv_nom_points + vv_total_syst_dev_high)*n_sf_vv

    vv_stat_high = vv_nom_points*(1+stat_uncertainty_rel)*n_sf_vv
    vv_stat_low = vv_nom_points*(1-stat_uncertainty_rel)*n_sf_vv

    vv_shape_high = vv_nom_points*(1+vv_shape_syst_points)*n_sf_vv
    vv_shape_low = vv_nom_points*(1-vv_shape_syst_points)*n_sf_vv

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.fill(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), color="y", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_shape_low, vv_shape_high), color="b", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color="r", alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_sf_vv, where="post", linestyle="dashed", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Diboson SF")

    ax.legend(['Fitted', 'Nominal', "Histogram Systematic", "Shape Systematic", "Statistical Systematic"])
    plt.savefig("plots/template_vv_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")

    ax.fill(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), color="y", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_shape_low, vv_shape_high), color="b", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color="r", alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_sf_vv, where="post", linestyle="--", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Diboson SF")
    ax.legend(['Fitted', 'Nominal', "Histogram Systematic", "Shape Systematic", "Statistical Systematic"])
    plt.savefig("plots/template_vv_sf_log.pdf")


    # W+Jets

    n_sf_wjets = ws.obj("n_sf_wjets").getVal()
    wjets_fitted = ws.obj("wjets_sf_sf_overallSyst_x_StatUncert_x_sf_wjets_syst_sf_ShapeSys")
    wjets_syst_nom = ws.obj("wjets_sf_sf_nominal")

    wjets_fitted_points = np.zeros(x.shape)
    wjets_nom_points = np.zeros(x.shape)
    wjets_shape_syst_points = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        wjets_fitted_points[i] = wjets_fitted.getVal()
        wjets_nom_points[i] = wjets_syst_nom.getVal()
        try:
            wjets_shape_syst_points[i] = ws.obj("gamma_wjets_syst_sf_bin_{0}_sigma".format(i)).getVal()
        except AttributeError:
            pass

    wjets_nom_points *= n_sf_wjets
    wjets_shape_syst_points *= wjets_nom_points

    wjets_syst_low = wjets_nom_points - wjets_shape_syst_points
    wjets_syst_high = wjets_nom_points + wjets_shape_syst_points

    wjets_stat_high = wjets_nom_points*(1+stat_uncertainty_rel)
    wjets_stat_low = wjets_nom_points*(1-stat_uncertainty_rel)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill(*get_step_fill_between(x, wjets_syst_low, wjets_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, wjets_stat_low, wjets_stat_high), color='r', alpha=0.5)
    ax.step(x, wjets_fitted_points, color="k", where="post")
    ax.step(x, wjets_nom_points, color="b", linestyle="--", where="post")
    ax.set_title("Non-prompt SF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template_wjets_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill(*get_step_fill_between(x, wjets_syst_low, wjets_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, wjets_stat_low, wjets_stat_high), color='r', alpha=0.5)
    ax.step(x, wjets_fitted_points, color="k", where="post")
    ax.step(x, wjets_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 500)
    ax.set_title("Non-prompt SF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template_wjets_sf_log.pdf")


    # Z
    n_sf_z = ws.obj("n_sf_z").getVal()
    z_fitted = ws.obj("z_sf_sf_overallSyst_x_StatUncert_x_sf_z_syst_sf_ShapeSys")
    z_syst_nom = ws.obj("z_sf_sf_nominal")

    z_fitted_points = np.zeros(x.shape)
    z_nom_points = np.zeros(x.shape)
    z_shape_syst_points = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        z_fitted_points[i] = z_fitted.getVal()
        z_nom_points[i] = z_syst_nom.getVal()
        try:
            z_shape_syst_points[i] = ws.obj("gamma_z_syst_sf_bin_{0}_sigma".format(i)).getVal()
        except AttributeError:
            pass

    z_nom_points *= n_sf_z
    z_shape_syst_points *= z_nom_points

    z_syst_low = z_nom_points - z_shape_syst_points
    z_syst_high = z_nom_points + z_shape_syst_points

    z_stat_high = z_nom_points*(1+stat_uncertainty_rel)
    z_stat_low = z_nom_points*(1-stat_uncertainty_rel)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill(*get_step_fill_between(x, z_syst_low, z_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, z_stat_low, z_stat_high), color='r', alpha=0.5)
    ax.step(x, z_fitted_points, color="k", where="post")
    ax.step(x, z_nom_points, color="b", linestyle="--", where="post")
    ax.set_title("Z SF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal',"Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template_z_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill(*get_step_fill_between(x, z_syst_low, z_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, z_stat_low, z_stat_high), color='r', alpha=0.5)
    ax.step(x, z_fitted_points, color="k", where="post")
    ax.step(x, z_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 1000)
    ax.set_title("Z SF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal',"Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template_z_sf_log.pdf")

def plot_fitted_of(ws):
    obs = ws.obj("obs_x_of")
    x = np.arange(10, 300, 10)


    # statistical uncertainties
    stat_uncertainty_means = np.zeros(x.shape)
    for i in xrange(len(x)):
        try:
            stat_uncertainty_means[i] = ws.obj("gamma_stat_of_bin_{0}_tau".format(i)).getVal()
        except AttributeError:
            pass

    stat_uncertainty_rel = 1./np.sqrt(stat_uncertainty_means)

    # top
    par = R.RooArgList(ws.obj("n_top_sf"), ws.obj("alpha_top_ratio"))
    scale = ws.obj("n_top_of_scale").getVal()
    n_top_of = R.RooFormulaVar("n_top_of", "n_top_of", "n_top_sf*{0}*(1+0.1*alpha_top_ratio)".format(scale), par)
    n_of_top = n_top_of.getVal()
    top_fitted = ws.obj("top_of_of_overallSyst_x_StatUncert")
    top_nominal = ws.obj("top_of_of_nominal")

    # import IPython
    # IPython.embed()

    top_fitted_points = np.zeros(x.shape)
    top_nom_points = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        top_fitted_points[i] = top_fitted.getVal()
        top_nom_points[i] = top_nominal.getVal()

    top_nom_points *= n_of_top
    top_stat_high = top_nom_points*(1+stat_uncertainty_rel)
    top_stat_low = top_nom_points*(1-stat_uncertainty_rel)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill(*get_step_fill_between(x, top_stat_low, top_stat_high), color='r', alpha=0.5)
    ax.step(x, top_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, top_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Top OF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal',"Statistical Systematic"])

    plt.savefig("plots/template_top_of.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill(*get_step_fill_between(x, top_stat_low, top_stat_high), color='r', alpha=0.5)
    ax.step(x, top_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, top_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_ylim(0.001, 5000)
    ax.set_title("Top OF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Statistical Systematic"])

    plt.savefig("plots/template_top_of_log.pdf")

    # VV

    # fitted shape
    # includes both fitted statistical and systematics
    # might need to multiply by bin width
    par = R.RooArgList(ws.obj("n_vv_sf"), ws.obj("alpha_vv_ratio"))
    scale = ws.obj("n_vv_of_scale").getVal()
    n_vv_of = R.RooFormulaVar("n_vv_of", "n_vv_of", "n_vv_sf*{0}*(1+0.1*alpha_vv_ratio)".format(scale), par)
    n_of_vv = n_vv_of.getVal()
    vv_fitted = ws.obj("vv_of_of_overallSyst_x_StatUncert_x_of_ww_syst_of_ShapeSys")
    vv_syst_nom = ws.obj("vv_of_of_Hist_alphanominal")

    n_vv_syst = 4
    vv_systs_low = []
    vv_systs_high = []

    for i in xrange(n_vv_syst):
        vv_systs_low.append(ws.obj("vv_of_of_Hist_alpha_{0}low".format(i)))
        vv_systs_high.append(ws.obj("vv_of_of_Hist_alpha_{0}high".format(i)))


    vv_fitted_points = np.zeros(x.shape)
    vv_nom_points = np.zeros(x.shape)

    vv_shape_syst_points = np.zeros(x.shape)

    vv_syst_low_points = np.zeros((n_vv_syst, len(x)))
    vv_syst_high_points = np.zeros((n_vv_syst, len(x)))

    for i in xrange(len(x)):
        obs.setVal(x[i])
        vv_fitted_points[i] = vv_fitted.getVal()
        vv_nom_points[i] = vv_syst_nom.getVal()
        try:
            vv_shape_syst_points[i] = ws.obj("gamma_ww_syst_of_bin_{0}_sigma".format(i)).getVal()
        except AttributeError:
            pass
        for j in xrange(n_vv_syst):
            vv_syst_low_points[j,i] = vv_systs_low[j].getVal()
            vv_syst_high_points[j,i] = vv_systs_high[j].getVal()


    # now get the deviations from the nominal
    vv_syst_dev_low = vv_syst_low_points - vv_nom_points
    vv_syst_dev_high = vv_syst_high_points - vv_nom_points

    vv_total_syst_dev_low = np.sqrt(np.sum(vv_syst_dev_low**2, axis=0))
    vv_total_syst_dev_high = np.sqrt(np.sum(vv_syst_dev_high**2, axis=0))


    vv_syst_low_points =  (vv_nom_points - vv_total_syst_dev_low)*n_of_vv
    vv_syst_high_points = (vv_nom_points + vv_total_syst_dev_high)*n_of_vv

    vv_stat_high = vv_nom_points*(1+stat_uncertainty_rel)*n_of_vv
    vv_stat_low = vv_nom_points*(1-stat_uncertainty_rel)*n_of_vv

    vv_shape_high = vv_nom_points*(1+vv_shape_syst_points)*n_of_vv
    vv_shape_low = vv_nom_points*(1-vv_shape_syst_points)*n_of_vv

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.fill(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), color="y", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_shape_low, vv_shape_high), color="b", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color="r", alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_of_vv, where="post", linestyle="dashed", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Diboson OF")

    ax.legend(['Fitted', 'Nominal', "Histogram Systematic", "Shape Systematic", "Statistical Systematic"])
    plt.savefig("plots/template_vv_of.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")

    ax.fill(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), color="y", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_shape_low, vv_shape_high), color="b", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color="r", alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_of_vv, where="post", linestyle="--", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Diboson OF")
    ax.legend(['Fitted', 'Nominal', "Histogram Systematic", "Shape Systematic", "Statistical Systematic"])
    plt.savefig("plots/template_vv_of_log.pdf")


    # W+Jets

    n_of_wjets = ws.obj("n_of_wjets").getVal()
    wjets_fitted = ws.obj("wjets_of_of_overallSyst_x_StatUncert_x_of_wjets_syst_of_ShapeSys")
    wjets_syst_nom = ws.obj("wjets_of_of_nominal")

    wjets_fitted_points = np.zeros(x.shape)
    wjets_nom_points = np.zeros(x.shape)
    wjets_shape_syst_points = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        wjets_fitted_points[i] = wjets_fitted.getVal()
        wjets_nom_points[i] = wjets_syst_nom.getVal()
        try:
            wjets_shape_syst_points[i] = ws.obj("gamma_wjets_syst_of_bin_{0}_sigma".format(i)).getVal()
        except AttributeError:
            pass

    wjets_nom_points *= n_of_wjets
    wjets_shape_syst_points *= wjets_nom_points

    wjets_syst_low = wjets_nom_points - wjets_shape_syst_points
    wjets_syst_high = wjets_nom_points + wjets_shape_syst_points

    wjets_stat_high = wjets_nom_points*(1+stat_uncertainty_rel)
    wjets_stat_low = wjets_nom_points*(1-stat_uncertainty_rel)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill(*get_step_fill_between(x, wjets_syst_low, wjets_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, wjets_stat_low, wjets_stat_high), color='r', alpha=0.5)
    ax.step(x, wjets_fitted_points, color="k", where="post")
    ax.step(x, wjets_nom_points, color="b", linestyle="--", where="post")
    ax.set_title("Non-prompt OF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template_wjets_of.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill(*get_step_fill_between(x, wjets_syst_low, wjets_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, wjets_stat_low, wjets_stat_high), color='r', alpha=0.5)
    ax.step(x, wjets_fitted_points, color="k", where="post")
    ax.step(x, wjets_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 500)
    ax.set_title("Non-prompt OF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template_wjets_of_log.pdf")


    # Z
    n_of_z = ws.obj("n_of_z").getVal()
    z_fitted = ws.obj("z_of_of_overallSyst_x_StatUncert_x_of_z_syst_of_ShapeSys")
    z_syst_nom = ws.obj("z_of_of_nominal")

    z_fitted_points = np.zeros(x.shape)
    z_nom_points = np.zeros(x.shape)
    z_shape_syst_points = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        z_fitted_points[i] = z_fitted.getVal()
        z_nom_points[i] = z_syst_nom.getVal()
        try:
            z_shape_syst_points[i] = ws.obj("gamma_z_syst_of_bin_{0}_sigma".format(i)).getVal()
        except AttributeError:
            pass

    z_nom_points *= n_of_z
    z_shape_syst_points *= z_nom_points

    z_syst_low = z_nom_points - z_shape_syst_points
    z_syst_high = z_nom_points + z_shape_syst_points

    z_stat_high = z_nom_points*(1+stat_uncertainty_rel)
    z_stat_low = z_nom_points*(1-stat_uncertainty_rel)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill(*get_step_fill_between(x, z_syst_low, z_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, z_stat_low, z_stat_high), color='r', alpha=0.5)
    ax.step(x, z_fitted_points, color="k", where="post")
    ax.step(x, z_nom_points, color="b", linestyle="--", where="post")
    ax.set_title("Z OF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template_z_of.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill(*get_step_fill_between(x, z_syst_low, z_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, z_stat_low, z_stat_high), color='r', alpha=0.5)
    ax.step(x, z_fitted_points, color="k", where="post")
    ax.step(x, z_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 1000)
    ax.set_title("Z OF")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template_z_of_log.pdf")

if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)
    print args
    sys.argv = [sys.argv[0], "-b"]


    file_name = args['<filename>']
    out_file = args['<outfile>']

    cut_val = float(args['--cut'])
    if cut_val < 0:
        cut = None
    else:
        cut = cut_val

    ncpu = int(args['--ncpu'])

    get_p = bool(args['-p'])
    do_minos = bool(args['-m'])

    res = run_bonly_fit(file_name, out_file, ncpu, get_p, do_minos=do_minos, cut=cut)

