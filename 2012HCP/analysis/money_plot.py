import numpy as np
import json
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('text.latex', preamble = '\usepackage{amsmath}')

from matplotlib.pyplot import *
from matplotlib.font_manager import *
switch_backend("Agg")

from matplotlib.ticker import MultipleLocator

import sys
sys.path.append("../")
sys.path.append("../../")
from data import *
from CMSPyLibs.plot import *

bins=30
plotrange = (0,300)

def make_money_plot(flavor):
    bkgtpl = []
    bkgwtpl = []
    bkgltpl = []

    fontp = FontProperties(family="Helvetica", size=12)
    fontpb = FontProperties(family="Helvetica", size=12, weight="book")


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
    bkgltpl.append("Non-prompt")
    print "Fake Pred: ", sf*data[sd['wjets_mct_high'+flavor]].mctperp.count()

    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.01, 10000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')
    h = hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl, zorder=1, linewidth=0.5)
    he = hist_errorbars( data[sd['sig'+flavor]].mctperp, xerrs=False, bins=bins, range=plotrange)
    he.set_label("Data")

    # move data to top of legend
    handles, labels = fig.get_legend_handles_labels()
    handles.insert(0,handles.pop())
    labels.insert(0,labels.pop())

    legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=2)
    fig.set_axisbelow(False)

    minorticks = MultipleLocator(10)
    fig.xaxis.set_minor_locator(minorticks)

    all_bkg = Series()
    all_bkg_w = Series()
    for mct, w in zip(bkgtpl, bkgwtpl):
        all_bkg = all_bkg.append(mct)
        all_bkg_w = all_bkg_w.append(Series(w))


    fig2 = subplot2grid((4,1),(3,0))
    hist_ratio(data[sd['sig'+flavor]].mctperp, all_bkg, all_bkg_w, bins=bins, range=plotrange)
    axhline(1, color="k")
    fig2.set_ylim(0, 2)
    fig2.set_ylim(0,2)
    fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    figtext(0.12, 0.92, r"CMS Preliminary $\sqrt{\text{s}}=8\;\text{TeV},$\quad L$_{\text{int}}=9.2\;\text{fb}^{-1}$", color='k',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

    savefig("plots/money"+flavor+".pdf")

    return fig

