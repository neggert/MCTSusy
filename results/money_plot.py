"""
This file contains the function that is used to make the "money plot". It reads
background fit information from the json files, so make sure those are up-to-date.
"""


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
sys.path.append("import/")
sys.path.append("../import/")
from config.data import *
from CMSPyLibs.plot import *

bins=29
plotrange = (10,300)

channels = ['sf', 'of']

def make_money_plot():

    fontp = FontProperties(family="Helvetica", size=12)
    fontpb = FontProperties(family="Helvetica", size=12, weight="book")


    f = open("fit_results.json")
    results = json.load(f)
    f.close()

    figs = []
    for ch in channels:
        bkgtpl = []
        bkgwtpl = []
        bkgltpl = []

        bkgtpl.append( data[sd['top_ctrl_'+ch]].mctperp )
        data_norm = data[sd['top_mct_low_'+ch]].mctperp.count()
        est_events = float(results[ch]['top'][0])
        sf = est_events/data_norm
        bkgwtpl.append( sf*np.ones(data[sd['top_ctrl_'+ch]].mctperp.count()) )
        bkgltpl.append("Top")

        bkgtpl.append( mc[smc['sig_'+ch]&(mc.mctype.isin(['WWW','WWG','WWTo2L2Nu']))].mctperp )
        allvv_norm = sum(mc[smc['sig_mct_low_'+ch]&((mc.mc_cat=='WV') | (mc.mc_cat=='ZZ'))].weight)
        est_events = float(results[ch]['vv'][0])
        sf = est_events/allvv_norm
        bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mctype.isin(['WWW','WWG','WWTo2L2Nu']))].weight)
        bkgltpl.append("WW")

        bkgtpl.append( mc[smc['sig_'+ch]&(mc.mctype=='WZTo3LNu')].mctperp )
        allvv_norm = sum(mc[smc['sig_mct_low_'+ch]&((mc.mc_cat=='WV') | (mc.mc_cat=='ZZ'))].weight)
        est_events = float(results[ch]['vv'][0])
        sf = est_events/allvv_norm
        bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mctype=='WZTo3LNu')].weight)
        bkgltpl.append("WZ")
        print "WZ Pred: ", sf*sum(mc[smc['sig_mct_high_'+ch]&(mc.mctype=='WZTo3LNu')].weight)

        bkgtpl.append( mc[smc['sig_'+ch]&(mc.mctype=='ZZTo2L2Nu')].mctperp )
        allvv_norm = sum(mc[smc['sig_mct_low_'+ch]&((mc.mc_cat=='WV') | (mc.mc_cat=='ZZ'))].weight)
        est_events = float(results[ch]['vv'][0])
        sf = est_events/allvv_norm
        bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mctype=='ZZTo2L2Nu')].weight)
        bkgltpl.append("ZZ")
        print "ZZ Pred: ", sf*sum(mc[smc['sig_mct_high_'+ch]&(mc.mctype=='ZZTo2L2Nu')].weight)

        bkgtpl.append( mc[smc['sig_'+ch]&(mc.mc_cat=="DY")].mctperp )
        data_norm = sum(mc[smc['sig_mct_low_'+ch]&(mc.mc_cat=="DY")].weight)
        est_events = float(results[ch]['z'][0])
        sf = est_events/data_norm
        bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mc_cat=="DY")].weight )
        bkgltpl.append("Z/$\gamma^*$")

        bkgtpl.append( data[sd['wjets_ctrl_'+ch]].mctperp )
        data_norm = data[sd['wjets_mct_low_'+ch]].mctperp.count()
        est_events = float(results[ch]['wjets'][0])
        sf = est_events/data_norm
        bkgwtpl.append( sf*np.ones(data[sd['wjets_ctrl_'+ch]].mctperp.count()) )
        bkgltpl.append("Non-prompt")

        f = figure(figsize=(6,6))
        f.set_facecolor('w')
        fig = subplot2grid((4,1),(0,0), rowspan=3)
        fig.set_yscale('log', nonposy='clip')
        fig.set_ylim(0.001, 1000)
        fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')
        h = hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl, zorder=1, linewidth=0.5)
        he = hist_errorbars( data[sd['sig_'+ch]].mctperp, xerrs=False, bins=bins, range=plotrange)
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
        hist_ratio(data[sd['sig_'+ch]].mctperp, all_bkg, all_bkg_w, bins=bins, range=plotrange)
        axhline(1, color="k")
        fig2.set_ylim(0.5,1.5)
        fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

        xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

        figtext(0.12, 0.92, r"CMS Preliminary $\sqrt{\text{s}}=8\;\text{TeV},$\quad L$_{\text{int}}=9.2\;\text{fb}^{-1}$", color='k',
                 fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

        savefig("plots/money"+ch+".pdf")

        figs.append(fig)

    return fig

if __name__ == '__main__':
    make_money_plot()
