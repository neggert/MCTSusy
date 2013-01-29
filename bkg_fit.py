#!/usr/bin/env python

#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    bkg_fit.py <signal_file> <mass1> <mass2> [-p] [--ncpu=<c>] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    --ncpu=<c>       Number of CPUs to use [default: 1]
    --channels=<c1,c2> Channels to use [default: of,sf]
    -p               Get p value (slow)

"""

from set_limits import *
from collections import defaultdict
import json
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import IndexLocator, FixedFormatter
import numpy as np

bkgs = ['top', 'vv', 'wjets', 'z']

def run_bonly_fit(sig_file, mass1, mass2, chans, ncpu, get_p, data_prefix="data", data_file_name="data.root"):
    prefix = "limits/"+sig_file[:-5]+"_{}_{}".format(mass1, mass2)

    try:
        create_histfactory(sig_file, prefix, mass1, mass2, chans, data_prefix=data_prefix, data_file_name=data_file_name)
    except:
        return None

    rfile = R.TFile(prefix+"_combined_meas_model.root")

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    model.GetParametersOfInterest().first().setVal(0.)
    model.GetParametersOfInterest().first().setConstant()

    pars = model.GetNuisanceParameters()

    err_pars = R.RooArgSet(pars.find("n_of_top"), pars.find("n_of_vv"), pars.find("n_of_z"), pars.find("n_of_wjets"),
                       pars.find("n_sf_top"), pars.find("n_sf_vv"), pars.find("n_sf_z"), pars.find("n_sf_wjets"))

    # run the fit
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR)
    res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save())

    fitPars = res.floatParsFinal()

    fitresults = defaultdict(dict)
    for ch in chans:
        for b in bkgs:
            fitvar = fitPars.find('n_{}_{}'.format(ch, b))
            fitresults[ch][b] = (fitvar.getVal(), fitvar.getError())

    f = open("fit_results.json", 'w')

    json.dump(fitresults, f, indent=3)

    f.close()

    # plot the relevant portion of the correlation matrix
    cor = res.correlationMatrix()
    cor = cor.GetSub(97, 104, 97, 104)
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


    for i in xrange(cor.GetNrows()):
        for j in xrange(cor.GetNcols()):
            c = cor[i][j]
            if c > 0: color='white'
            else: color='black'
            size = np.sqrt(np.abs(c))
            rect = Rectangle([i-size/2, j-size/2], size, size, facecolor=color, edgecolor='black')
            ax.add_patch(rect)
    ax.autoscale_view()
    plt.gca().invert_yaxis()
    plt.tight_layout()

    plt.savefig("plots/correlation.pdf")
    raw_input("...")

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



if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)

    m1 = int(args['<mass1>'])
    m2 = int(args['<mass2>'])

    sig_file = args['<signal_file>']

    try:
        chans = args['--channels'].split(",")
    except AttributeError:
        chans = ['of', 'sf']

    ncpu = int(args['--ncpu'])

    get_p = bool(args['-p'])

    res = run_bonly_fit(sig_file, m1, m2, chans, ncpu, get_p)

