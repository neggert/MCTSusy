import sys
sys.path.append("../../")
from CMSPyLibs.plot import *
from data import *

from matplotlib.pyplot import *

def plot_fake_closure(flavor):
    faketruth = smc['sig'+flavor] & (mc.mc_cat=='fake')

    nbins=5
    prange=(1,101)

    figure(figsize=(14,7))
    fig = subplot(121)
    fig.set_yscale('log', nonposy='clip')
    hist( mc[faketruth].mctperp, weights=mc[faketruth].weight, bins=nbins, range=prange, histtype="step", stacked=True,\
        normed=True, label="W+Jets MC Truth", color='#440088')
    he = hist_errorbars( mc[smc['wjets_ctrl'+flavor]].mctperp.values, weights=mc[smc['wjets_ctrl'+flavor]].weight.values, bins=nbins, range=prange, normed=True,\
        xerrs=False, label="Control", color='k')
    he.set_label("Control Region MC")
    ylim(1.e-4, .1)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    subplot(122)
    hist( mc[faketruth].mctperp, weights=mc[faketruth].weight, bins=nbins, range=prange, histtype="step", stacked=True,\
        normed=True, label="W+Jets MC Truth", color='#440088')
    he = hist_errorbars( mc[smc['wjets_ctrl'+flavor]].mctperp.values, weights=mc[smc['wjets_ctrl'+flavor]].weight.values, bins=nbins, range=prange, normed=True,\
        xerrs=False, label="Control Region MC", color='k')
    he.set_label("Control Region MC")
    ylim(0, 0.04)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")

    savefig("plots/closure_fake{}.pdf".format(flavor))

