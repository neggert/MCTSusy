#!/usr/bin/env python

#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    bkg_fit.py <filename> [-pm] [--ncpu=<c>] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    --ncpu=<c>       Number of CPUs to use [default: 1]
    --channels=<c1,c2> Channels to use [default: of,of]
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

bkgs = ['top', 'vv', 'wjets', 'z']

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

    return fill_x, fill_y1, fill_y2

def run_bonly_fit(file_name, ncpu, get_p, data_prefix="data", data_file_name="data.root", do_minos=False):

    rfile = R.TFile(file_name)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    model.GetParametersOfInterest().first().setVal(0.)
    model.GetParametersOfInterest().first().setConstant()

    pars = model.GetNuisanceParameters()

    err_pars = R.RooArgSet(pars.find("n_of_top"), pars.find("n_of_vv"), pars.find("n_of_z"), pars.find("n_of_wjets"),
                       pars.find("n_of_top"), pars.find("n_of_vv"), pars.find("n_of_z"), pars.find("n_of_wjets"))

    # run the fit
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.PrintLevel(0))
    if do_minos:
        res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.Minos())

    fitPars = res.floatParsFinal()

    nPar = fitPars.getSize()

    t = PrettyTable()
    if do_minos:
        t.field_names = ['', "Value", "Parabolic Error", "Minos Down", "Minos Up"]
    else:
        t.field_names = ['', "Value", "Parabolic Error"]
    t.vertical_char = "&"
    for i in xrange(nPar):
        name = fitPars.at(i).GetName()
        val =  fitPars.at(i).getVal()
        err =  fitPars.at(i).getError()
        if do_minos:
            minos_up = fitPars.at(i).getErrorLo()
            minos_down = fitPars.at(i).getErrorHi()
            t.add_row((name, "{0:.3g}".format(val), "{0:.3g}".format(err), "{0:.3g}".format(minos_down), "{0:.3g}".format(minos_up)))

        else:
            t.add_row((name, "{0:.3g}".format(val), "{0:.3g}".format(err)))

    print t

    fitresults = defaultdict(dict)
    chans = ['of','sf']
    for ch in chans:
        for b in bkgs:
            fitvar = fitPars.find('n_{}_{}'.format(ch, b))
            if b == 'z':
                b = "DY"
            elif b == 'wjets':
                b = 'fake'
            fitresults[ch][b] = (fitvar.getVal(), fitvar.getError())

    f = open("fit_results.json", 'w')

    json.dump(fitresults, f, indent=3)

    f.close()

    # plot the relevant portion of the correlation matrix
    fullcor = res.correlationMatrix()
    cor = fullcor.GetSub(96, 103, 96, 103)
    # import pdb; pdb.set_trace()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.patch.set_facecolor('gray')
    ax.set_aspect('equal', 'box')
    labels = FixedFormatter(['',
                            'OF Top',
                             'OF Diboson',
                             'OF Z',
                             'OF WJets',
                             'SF Top',
                             'SF Diboson',
                             'SF Z',
                             'SF WJets'])
    ax.xaxis.set_major_formatter(labels)
    ax.xaxis.tick_top()
    ax.yaxis.set_major_formatter(labels)
    for t in ax.get_xticklabels():
        t.set_rotation(40)
        t.set_ha('left')


    # make Hinton-style correlation plot
    for i in xrange(cor.GetNrows()):
        for j in xrange(cor.GetNcols()):
            # if i<=j: continue
            c = cor[i][j]
            if abs(c) < 0.01: continue
            if c > 0: color='white'
            else: color='black'
            size = np.sqrt(np.abs(c))
            rect = Rectangle([i-size/2, j-size/2], size, size, facecolor=color, edgecolor='black')
            ax.add_patch(rect)
    ax.autoscale_view()
    plt.gca().invert_yaxis()
    plt.tight_layout()

    plt.savefig("plots/correlation.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.patch.set_facecolor('gray')
    ax.set_aspect('equal', 'box')

    cor.Print()


    # make Hinton-style correlation plot
    for i in xrange(fullcor.GetNrows()):
        for j in xrange(fullcor.GetNcols()):
            # if i<=j: continue
            c = fullcor[i][j]
            if abs(c) < 0.01: continue
            if c > 0: color='white'
            else: color='black'
            size = np.sqrt(np.abs(c))
            rect = Rectangle([i-size/2, j-size/2], size, size, facecolor=color, edgecolor='black')
            ax.add_patch(rect)
    ax.set_xlim(0, 140)
    ax.set_ylim(0, 140)
    plt.gca().invert_yaxis()
    plt.tight_layout()

    plt.savefig("plots/correlation_full.pdf")
    # raw_input("...")

    model.SetSnapshot(model.GetParametersOfInterest())

    # plot_fitted_sf(ws)
    # plot_fitted_of(ws)


    raw_input("...")


    # plot the fitted templates

    # get the test statistic on data    
    R.gROOT.ProcessLineSync(".L KS/AndersonDarlingTestStat.cc+")
    AD = R.RooStats.AndersonDarlingTestStat(ws.obj("simPdf"), model.GetPdf())
    ts = AD.Evaluate(data, model.GetParametersOfInterest())

    # import IPython
    # IPython.embed()

    # calculate a p-value
    if get_p:

        sampler = R.RooStats.ToyMCSampler(AD, 10)
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
    n_sf_top = ws.obj("n_sf_top").getVal()
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
    ax.fill_between(*get_step_fill_between(x, top_stat_low, top_stat_high), color='r', alpha=0.5)
    ax.step(x, top_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, top_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Top SF")
    plt.savefig("plots/template_top_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill_between(*get_step_fill_between(x, top_stat_low, top_stat_high), color='r', alpha=0.5)
    ax.step(x, top_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, top_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_ylim(0.001, 5000)
    ax.set_title("Top SF")
    plt.savefig("plots/template_top_sf_log.pdf")

    # VV

    # fitted shape
    # includes both fitted statistical and systematics
    # might need to multiply by bin width
    n_sf_vv = ws.obj("n_sf_vv").getVal()
    vv_fitted = ws.obj("vv_sf_sf_overallSyst_x_StatUncert")
    vv_syst_nom = ws.obj("vv_sf_sf_Hist_alphanominal")

    n_vv_syst = 4
    vv_systs_low = []
    vv_systs_high = []

    for i in xrange(n_vv_syst):
        vv_systs_low.append(ws.obj("vv_sf_sf_Hist_alpha_{0}low".format(i)))
        vv_systs_high.append(ws.obj("vv_sf_sf_Hist_alpha_{0}high".format(i)))


    vv_fitted_points = np.zeros(x.shape)
    vv_nom_points = np.zeros(x.shape)

    vv_syst_low_points = np.zeros((n_vv_syst, len(x)))
    vv_syst_high_points = np.zeros((n_vv_syst, len(x)))

    for i in xrange(len(x)):
        obs.setVal(x[i])
        vv_fitted_points[i] = vv_fitted.getVal()
        vv_nom_points[i] = vv_syst_nom.getVal()
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

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.fill_between(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color='r', alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_sf_vv, where="post", linestyle="dashed", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Diboson SF")
    plt.savefig("plots/template_vv_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")

    ax.fill_between(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color='r', alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_sf_vv, where="post", linestyle="--", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Diboson SF")
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
    ax.fill_between(*get_step_fill_between(x, wjets_syst_low, wjets_syst_high), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, wjets_stat_low, wjets_stat_high), color='r', alpha=0.5)
    ax.step(x, wjets_fitted_points, color="k", where="post")
    ax.step(x, wjets_nom_points, color="b", linestyle="--", where="post")
    ax.set_title("Non-prompt SF")
    plt.savefig("plots/template_wjets_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill_between(*get_step_fill_between(x, wjets_syst_low, wjets_syst_high), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, wjets_stat_low, wjets_stat_high), color='r', alpha=0.5)
    ax.step(x, wjets_fitted_points, color="k", where="post")
    ax.step(x, wjets_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 500)
    ax.set_title("Non-prompt SF")
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
    ax.fill_between(*get_step_fill_between(x, z_syst_low, z_syst_high), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, z_stat_low, z_stat_high), color='r', alpha=0.5)
    ax.step(x, z_fitted_points, color="k", where="post")
    ax.step(x, z_nom_points, color="b", linestyle="--", where="post")
    ax.set_title("Z SF")
    plt.savefig("plots/template_z_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill_between(*get_step_fill_between(x, z_syst_low, z_syst_high), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, z_stat_low, z_stat_high), color='r', alpha=0.5)
    ax.step(x, z_fitted_points, color="k", where="post")
    ax.step(x, z_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 1000)
    ax.set_title("Z SF")
    plt.savefig("plots/template_z_sf_log.pdf")

    plt.show()

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
    n_of_top = ws.obj("n_of_top").getVal()
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
    ax.fill_between(*get_step_fill_between(x, top_stat_low, top_stat_high), color='r', alpha=0.5)
    ax.step(x, top_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, top_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Top OF")
    plt.savefig("plots/template_top_of.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill_between(*get_step_fill_between(x, top_stat_low, top_stat_high), color='r', alpha=0.5)
    ax.step(x, top_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, top_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_ylim(0.001, 5000)
    ax.set_title("Top OF")
    plt.savefig("plots/template_top_of_log.pdf")

    # VV

    # fitted shape
    # includes both fitted statistical and systematics
    # might need to multiply by bin width
    n_of_vv = ws.obj("n_of_vv").getVal()
    vv_fitted = ws.obj("vv_of_of_overallSyst_x_StatUncert")
    vv_syst_nom = ws.obj("vv_of_of_Hist_alphanominal")

    n_vv_syst = 4
    vv_systs_low = []
    vv_systs_high = []

    for i in xrange(n_vv_syst):
        vv_systs_low.append(ws.obj("vv_of_of_Hist_alpha_{0}low".format(i)))
        vv_systs_high.append(ws.obj("vv_of_of_Hist_alpha_{0}high".format(i)))


    vv_fitted_points = np.zeros(x.shape)
    vv_nom_points = np.zeros(x.shape)

    vv_syst_low_points = np.zeros((n_vv_syst, len(x)))
    vv_syst_high_points = np.zeros((n_vv_syst, len(x)))

    for i in xrange(len(x)):
        obs.setVal(x[i])
        vv_fitted_points[i] = vv_fitted.getVal()
        vv_nom_points[i] = vv_syst_nom.getVal()
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

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.fill_between(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color='r', alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_of_vv, where="post", linestyle="dashed", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Diboson OF")
    plt.savefig("plots/template_vv_of.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")

    ax.fill_between(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color='r', alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_of_vv, where="post", linestyle="--", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_ylim(0.001, 5000)
    ax.set_title("Diboson OF")
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
    ax.fill_between(*get_step_fill_between(x, wjets_syst_low, wjets_syst_high), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, wjets_stat_low, wjets_stat_high), color='r', alpha=0.5)
    ax.step(x, wjets_fitted_points, color="k", where="post")
    ax.step(x, wjets_nom_points, color="b", linestyle="--", where="post")
    ax.set_title("Non-prompt OF")
    plt.savefig("plots/template_wjets_of.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill_between(*get_step_fill_between(x, wjets_syst_low, wjets_syst_high), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, wjets_stat_low, wjets_stat_high), color='r', alpha=0.5)
    ax.step(x, wjets_fitted_points, color="k", where="post")
    ax.step(x, wjets_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 500)
    ax.set_title("Non-prompt OF")
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
    ax.fill_between(*get_step_fill_between(x, z_syst_low, z_syst_high), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, z_stat_low, z_stat_high), color='r', alpha=0.5)
    ax.step(x, z_fitted_points, color="k", where="post")
    ax.step(x, z_nom_points, color="b", linestyle="--", where="post")
    ax.set_title("Z OF")
    plt.savefig("plots/template_z_of.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill_between(*get_step_fill_between(x, z_syst_low, z_syst_high), alpha=0.5)
    ax.fill_between(*get_step_fill_between(x, z_stat_low, z_stat_high), color='r', alpha=0.5)
    ax.step(x, z_fitted_points, color="k", where="post")
    ax.step(x, z_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 500)
    ax.set_title("Z OF")
    plt.savefig("plots/template_z_of_log.pdf")

    plt.show()

if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)
    print args

    file_name = args['<filename>']



    ncpu = int(args['--ncpu'])

    get_p = bool(args['-p'])
    do_minos = bool(args['-m'])

    res = run_bonly_fit(file_name, ncpu, get_p, do_minos=do_minos)

