import numpy as np
import json
from matplotlib.pyplot import *

import sys
sys.path.append("../")
sys.path.append("../../")
from data import *
from CMSPyLibs.plot import *

import ROOT

bins=30
plotrange = (0,300)

def make_money_plot_root(flavor):
    bkgtpl = []
    bkgwtpl = []
    bkgltpl = []

    f = open("results"+flavor+".json")
    results = json.load(f)
    f.close()

    bkgtpl.append( data[sd['top_ctrl'+flavor]].mctperp )
    data_norm = data[sd['top_mct_low'+flavor]].mctperp.count()
    est_events = float(results['low']['Top'][0])
    sf = est_events/data_norm
    bkgwtpl.append( sf*np.ones(data[sd['top_ctrl'+flavor]].mctperp.count()) )
    bkgltpl.append("Top")
    print "Top Pred: ", sf*data[sd['top_mct_high'+flavor]].mctperp.count()

    bkgtpl.append( mc[smc['sig'+flavor]&(mc.mctype.isin(['WWW','WWG','WWTo2L2Nu']))].mctperp )
    data_norm = sum(mc[smc['sig_mct_low'+flavor]&(mc.mctype.isin(['WWW','WWG','WWTo2L2Nu']))].weight)
    est_events = float(results['low']['WW'][0])
    sf = est_events/data_norm
    bkgwtpl.append( sf*mc[smc['sig'+flavor]&(mc.mctype.isin(['WWW','WWG','WWTo2L2Nu']))].weight)
    bkgltpl.append("WW")
    print "WW Pred: ", sf*sum(mc[smc['sig_mct_high'+flavor]&(mc.mctype.isin(['WWW','WWG','WWTo2L2Nu']))].weight)

    bkgtpl.append( mc[smc['sig'+flavor]&(mc.mctype=='WZTo3LNu')].mctperp )
    data_norm = sum(mc[smc['sig_mct_low'+flavor]&(mc.mctype=='WZTo3LNu')].weight)
    est_events = float(results['low']['WZ'][0])
    sf = est_events/data_norm
    bkgwtpl.append( sf*mc[smc['sig'+flavor]&(mc.mctype=='WZTo3LNu')].weight)
    bkgltpl.append("WZ")
    print "WZ Pred: ", sf*sum(mc[smc['sig_mct_high'+flavor]&(mc.mctype=='WZTo3LNu')].weight)

    bkgtpl.append( mc[smc['sig'+flavor]&(mc.mctype=='ZZTo2L2Nu')].mctperp )
    data_norm = sum(mc[smc['sig_mct_low'+flavor]&(mc.mctype=='ZZTo2L2Nu')].weight)
    est_events = float(results['low']['ZZ'][0])
    sf = est_events/data_norm
    bkgwtpl.append( sf*mc[smc['sig'+flavor]&(mc.mctype=='ZZTo2L2Nu')].weight)
    bkgltpl.append("ZZ")
    print "ZZ Pred: ", sf*sum(mc[smc['sig_mct_high'+flavor]&(mc.mctype=='ZZTo2L2Nu')].weight)

    bkgtpl.append( mc[smc['sig'+flavor]&(mc.mc_cat=="DY")].mctperp )
    data_norm = sum(mc[smc['sig_mct_low'+flavor]&(mc.mc_cat=="DY")].weight)
    est_events = float(results['low']['DY'][0])
    sf = est_events/data_norm
    bkgwtpl.append( sf*mc[smc['sig'+flavor]&(mc.mc_cat=="DY")].weight )
    bkgltpl.append("Z/$\gamma^*$")
    print "DY Pred: ", sf*sum(mc[smc['sig_mct_high'+flavor]&(mc.mc_cat=="DY")].weight)

    bkgtpl.append( data[sd['wjets_ctrl'+flavor]].mctperp )
    data_norm = data[sd['wjets_mct_low'+flavor]].mctperp.count()
    est_events = float(results['low']['W'][0])
    sf = est_events/data_norm
    bkgwtpl.append( sf*np.ones(data[sd['wjets_ctrl'+flavor]].mctperp.count()) )
    bkgltpl.append("Fake")
    print "Fake Pred: ", sf*data[sd['wjets_mct_high'+flavor]].mctperp.count()


    # figure(figsize=(6,6))
    # fig = subplot2grid((4,1),(0,0), rowspan=3)
    # fig.set_yscale('log', nonposy='clip')
    # fig.set_ylim(0.01, 10000)
    # h = hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl, zorder=1)
    # he = hist_errorbars( data[sd['sig'+flavor]].mctperp, xerrs=False, bins=bins, range=plotrange)
    # he.set_label("Data")
    # legend()
    # fig.set_axisbelow(False)

    all_bkg = Series()
    all_bkg_w = Series()
    for mct, w in zip(bkgtpl, bkgwtpl):
        all_bkg = all_bkg.append(mct)
        all_bkg_w = all_bkg_w.append(Series(w))

    ROOT.gROOT.ProcessLineSync(".L tdrstyle.C")
    bkg_stack = ROOT.THStack("bkg", "bkg")
    colors = [ROOT.kBlue, ROOT.kGreen-1, ROOT.kRed, ROOT.kCyan, ROOT.kMagenta, ROOT.kYellow]
    i = 0
    l = ROOT.TLegend(.7, .5, .89, .7)
    hists = []
    for d, w, c in zip(bkgtpl, bkgwtpl, colors)[::-1]:
        i += 1
        h = ROOT.TH1D("hist"+str(i), "hist", bins, plotrange[0], plotrange[1])
        h.SetFillColor(c)
        h.SetLineColor(ROOT.kBlack)
        map(h.Fill, d, w)
        bkg_stack.Add(h)
        hists.append(h)

    c = ROOT.TCanvas("c1", "c1", 500, 500)
    c.cd()
    c.SetLogy()
    bkg_stack.SetMinimum(0.04)
    bkg_stack.SetMaximum(5000)
    bkg_stack.SetTitle("")
    ROOT.gStyle.SetOptStat(0)
    bkg_stack.Draw()
    # c.BuildLegend()

    names = ["Top", "WW", "WZ", "ZZ", "Z", "Fake"]
    for h, n in zip(hists[::-1], names):
        l.AddEntry(h, n, "F")

    l.SetFillColor(ROOT.kWhite)
    l.SetLineColor(ROOT.kWhite)
    l.Draw("same")

    text = ROOT.TLatex(120, 2000, "#splitline{CMS Preliminary}{#sqrt{s}=8 TeV, L_{int}=9.2 fb^{-1}}")
    text.Draw()


    raw_input("...")


    # fig2 = subplot2grid((4,1),(3,0), sharex=fig)
    # hist_ratio(data[sd['sig'+flavor]].mctperp, all_bkg, all_bkg_w, bins=bins, range=plotrange)
    # fig2.set_ylim(0.5, 1.5)
    # fig2.set_ylabel("Data/Prediction")

    # xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")

    # savefig("plots/money"+flavor+".pdf")


