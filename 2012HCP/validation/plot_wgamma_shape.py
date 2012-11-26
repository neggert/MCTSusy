import sys
sys.path.append("../../")
from CMSPyLibs.plot import *
from data import *

from matplotlib.pyplot import *

def plot_wgamma_shape(flavor):
    nbins=10
    prange=(0,200)

    wg_sig = mc[smc['sig'+flavor]&(mc.mc_cat=="wgamma")]
    wz_sig = mc[smc['sig'+flavor]&(mc.mctype=="WZTo3LNu")]

    figure(figsize=(14,7))
    fig = subplot(121)
    fig.set_yscale('log', nonposy='clip')
    hist( wz_sig.mctperp, weights=wz_sig.weight, bins=nbins, range=prange, histtype="step", stacked=True,\
        normed=True, label="WZ Shape", color='#440088')
    he = hist_errorbars( wg_sig.mctperp.values, weights=wg_sig.weight.values, bins=nbins, range=prange, normed=True,\
        xerrs=False, label="Control", color='k')
    he.set_label("WGstar shape")
    ylim(1.e-4, .1)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    subplot(122)
    hist( wz_sig.mctperp, weights=wz_sig.weight, bins=nbins, range=prange, histtype="step", stacked=True,\
        normed=True, label="WZ Shape", color='#440088')
    he = hist_errorbars( wg_sig.mctperp.values, weights=wg_sig.weight.values, bins=nbins, range=prange, normed=True,\
        xerrs=False, label="Control", color='k')
    he.set_label("WGstar shape")
    ylim(0, 0.04)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")

    savefig("plots/wgamm_shape{}.pdf".format(flavor))
