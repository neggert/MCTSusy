
import sys
sys.path.append("../../")
from CMSPyLibs.plot import *
from data import *

from matplotlib.pyplot import *

def plot_top_closure(flavor):
    """Make closure plot for the top control sample"""
    toptruth = smc['sig'+flavor] & (mc.mc_cat == 'top') &(mc.mctype!='tW')
    twtruth = smc['sig'+flavor] & (mc.mctype=='tW')

    tfig = figure(figsize=(14,7))
    nbins = 29
    nrange = (5., 150.)
    fig = subplot(121)
    fig.set_yscale('log', nonposy='clip')
    hist( mc[toptruth].mctperp, weights=mc[toptruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
        normed=True, label="Top MC Truth")
    hist( mc[twtruth].mctperp, weights=mc[twtruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
        normed=True, label="tW MC Truth")
    he = hist_errorbars( mc[smc['top_ctrl'+flavor]].mctperp.values, weights=mc[smc['top_ctrl'+flavor]].weight.values, bins=nbins, range=nrange, normed=True,\
        xerrs=False, color="k")
    he.set_label("$\geq$1 b-tags MC")
    ylim(1.e-7, .1)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    subplot(122)
    hist( mc[toptruth].mctperp, weights=mc[toptruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
        normed=True, label="Top MC Truth")
    hist( mc[twtruth].mctperp, weights=mc[twtruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
        normed=True, label="tW MC Truth")
    he = hist_errorbars( mc[smc['top_ctrl'+flavor]].mctperp.values, weights=mc[smc['top_ctrl'+flavor]].weight.values, bins=nbins, range=nrange, normed=True,\
        xerrs=False, color="k")
    he.set_label("$\geq$1 b-tags MC")
    ylim(0, 0.03)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    show()

    savefig("plots/closure_top{}.pdf".format(flavor))

