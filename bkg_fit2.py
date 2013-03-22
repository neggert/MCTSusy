#!/usr/bin/env python

#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    bkg_fit2.py <filename> [-pm] [--ncpu=<c>]

Options:
    -h --help        Show this screen.
    --ncpu=<c>       Number sf CPUs to use [default: 1]
    -p               Get p value (slow)
    -m               Run MINOS (slow)


"""


from set_limits2 import *
from collections import defaultdict
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import IndexLocator, FixedFormatter
from bkg_fit import get_step_fill_between
from prettytable import PrettyTable

bkgs = ['of', 'vv', 'wjets', 'z']

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

    # run the fit
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.PrintLevel(0))
    if do_minos:
        res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.Minos())

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

    fitresults = {}
    for b in bkgs:
        fitvar = fitPars.find('n_{}'.format(b))
        if b == 'z':
            b = "DY"
        elif b == 'wjets':
            b = 'fake'
        fitresults[b] = (fitvar.getVal(), fitvar.getError())

    f = open("fit_results2.json", 'w')

    json.dump(fitresults, f, indent=3)

    f.close()

    plot_fitted_sf(ws)

    fullcor = res.correlationMatrix()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.patch.set_facecolor('gray')
    ax.set_aspect('equal', 'box')

    # make Hinton-style correlation plot
    for i in xrange(fullcor.GetNrows()):
        for j in xrange(fullcor.GetNcols()):
            # if i<=j: continue
            c = fullcor[i][j]
            if abs(c) < 0.01: continue
            if c > 0: color='white'
            else: color='black'
            size = np.sqrt(np.abs(c))
            rect = Rectangle([i-size/2, j-size/2], size, size, facecolor=color, edgecolor='black', lw=0.1)
            ax.add_patch(rect)
    ax.set_xlim(0, fullcor.GetNrows())
    ax.set_ylim(0, fullcor.GetNrows())
    plt.gca().invert_yaxis()
    plt.tight_layout()

    plt.savefig("plots/correlation_full_sf.pdf")

    model.SetSnapshot(model.GetParametersOfInterest())

    R.gROOT.ProcessLineSync(".L KS/AndersonDarlingTestStat.cc+")

    AD = R.RooStats.AndersonDarlingTestStat(model.GetPdf())

    # get the test statistic on data
    ts = AD.Evaluate(data, model.GetParametersOfInterest())

    if get_p:

        sampler = R.RooStats.ToyMCSampler(AD, 500)
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
    n_fs = ws.obj("n_of")
    fs_fitted = ws.obj("of_sf_overallSyst_x_StatUncert")
    fs_nominal = ws.obj("of_sf_nominal")

    # import IPython
    # IPython.embed()

    fs_fitted_points = np.zeros(x.shape)
    fs_nom_points = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        fs_fitted_points[i] = fs_fitted.getVal()
        fs_nom_points[i] = fs_nominal.getVal()

    fs_nom_points *= n_fs.getVal()
    fs_stat_high = fs_nom_points*(1+stat_uncertainty_rel)
    fs_stat_low = fs_nom_points*(1-stat_uncertainty_rel)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill(*get_step_fill_between(x, fs_stat_low, fs_stat_high), color='r', alpha=0.5)
    ax.step(x, fs_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, fs_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Flavor-Symmetric")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal',"Statistical Systematic"])

    plt.savefig("plots/template2_fs.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill(*get_step_fill_between(x, fs_stat_low, fs_stat_high), color='r', alpha=0.5)
    ax.step(x, fs_nom_points, where="post", linestyle="dashed", color="b")
    ax.step(x, fs_fitted_points, where="post", color="k")
    ax.set_xlim(min(x), max(x)+10)
    ax.set_ylim(0.001, 5000)
    ax.set_title("Flavor-Symmetric")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Statistical Systematic"])

    plt.savefig("plots/template2_fs_log.pdf")

    # VV

    # fitted shape
    # includes both fitted statistical and systematics
    # might need to multiply by bin width
    n_vv_sf = ws.obj("n_vv")
    n_sf_vv = n_vv_sf.getVal()
    vv_fitted = ws.obj("vv_sf_overallSyst_x_StatUncert")
    vv_syst_nom = ws.obj("vv_sf_Hist_alphanominal")

    n_vv_syst = 2
    vv_systs_low = []
    vv_systs_high = []

    for i in xrange(n_vv_syst):
        vv_systs_low.append(ws.obj("vv_sf_Hist_alpha_{0}low".format(i)))
        vv_systs_high.append(ws.obj("vv_sf_Hist_alpha_{0}high".format(i)))


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

    ax.fill(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), color="y", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color="b", alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_sf_vv, where="post", linestyle="dashed", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Diboson")

    ax.legend(['Fitted', 'Nominal', "Histogram Systematic", "Statistical Systematic"])
    plt.savefig("plots/template2_vv_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")

    ax.fill(*get_step_fill_between(x, vv_syst_low_points, vv_syst_high_points), color="y", alpha=0.5)
    ax.fill(*get_step_fill_between(x, vv_stat_low, vv_stat_high), color="r", alpha=0.5)
    ax.step(x, vv_fitted_points, where="post", color="k")
    ax.step(x, vv_nom_points*n_sf_vv, where="post", linestyle="--", color="b")

    ax.set_xlim(min(x), max(x)+10)
    ax.set_title("Diboson")
    ax.legend(['Fitted', 'Nominal', "Histogram Systematic", "Statistical Systematic"])
    plt.savefig("plots/template2_vv_sf_log.pdf")


    # W+Jets

    n_sf_wjets = ws.obj("n_wjets").getVal()
    wjets_fitted = ws.obj("wjets_sf_overallSyst_x_StatUncert_x_sf_wjets_syst_ShapeSys")
    wjets_syst_nom = ws.obj("wjets_sf_nominal")

    wjets_fitted_points = np.zeros(x.shape)
    wjets_nom_points = np.zeros(x.shape)
    wjets_shape_syst_points = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        wjets_fitted_points[i] = wjets_fitted.getVal()
        wjets_nom_points[i] = wjets_syst_nom.getVal()
        try:
            wjets_shape_syst_points[i] = ws.obj("gamma_wjets_syst_bin_{0}_sigma".format(i)).getVal()
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
    ax.set_title("Non-prompt")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template2_wjets_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill(*get_step_fill_between(x, wjets_syst_low, wjets_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, wjets_stat_low, wjets_stat_high), color='r', alpha=0.5)
    ax.step(x, wjets_fitted_points, color="k", where="post")
    ax.step(x, wjets_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 500)
    ax.set_title("Non-prompt")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template2_wjets_sf_log.pdf")


    # Z
    n_sf_z = ws.obj("n_z").getVal()
    z_fitted = ws.obj("z_sf_overallSyst_x_StatUncert_x_sf_z_syst_ShapeSys")
    z_syst_nom = ws.obj("z_sf_nominal")

    z_fitted_points = np.zeros(x.shape)
    z_nom_points = np.zeros(x.shape)
    z_shape_syst_points = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        z_fitted_points[i] = z_fitted.getVal()
        z_nom_points[i] = z_syst_nom.getVal()
        try:
            z_shape_syst_points[i] = ws.obj("gamma_z_syst_bin_{0}_sigma".format(i)).getVal()
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
    ax.set_title("Z")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template2_z_sf.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log", nonposy="clip")
    ax.fill(*get_step_fill_between(x, z_syst_low, z_syst_high), alpha=0.5)
    ax.fill(*get_step_fill_between(x, z_stat_low, z_stat_high), color='r', alpha=0.5)
    ax.step(x, z_fitted_points, color="k", where="post")
    ax.step(x, z_nom_points, color="b", linestyle="--", where="post")
    ax.set_ylim(0.01, 1000)
    ax.set_title("Z")
    ax.set_xlabel(r"$M_{\mathrm{CT}\perp}$ (GeV)")
    ax.legend(['Fitted', 'Nominal', "Shape Systematic", "Statistical Systematic"])

    plt.savefig("plots/template2_z_sf_log.pdf")

if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)
    print args

    file_name = args['<filename>']



    ncpu = int(args['--ncpu'])

    get_p = bool(args['-p'])
    do_minos = bool(args['-m'])

    res = run_bonly_fit(file_name, ncpu, get_p, do_minos=do_minos)
