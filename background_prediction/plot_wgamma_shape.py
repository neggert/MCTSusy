import sys
sys.path.append("../../")
from CMSPyLibs.plot import *
from config.data import *

from matplotlib.pyplot import *

def plot_wgamma_shape(flavor):
    """Compare the shapes of the wgamma and WZ backgrounds"""
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

def plot_tt_ttw_shape(flavor):
    """Compare the shapes of the wgamma and WZ backgrounds"""
    nbins=29
    prange=(10,300)

    a_sig = mc[smc['sig'+flavor]&(mc.mctype=="ttbar")]
    b_sig = mc[smc['sig'+flavor]&(mc.mctype=="ttW")]

    figure(figsize=(14,7))
    fig = subplot(121)
    fig.set_yscale('log', nonposy='clip')
    hist( a_sig.mctperp, weights=a_sig.weight, bins=nbins, range=prange, histtype="step", stacked=True,\
        normed=True, label="Ttbar Shape", color='#440088')
    he = hist_errorbars( b_sig.mctperp.values, weights=b_sig.weight.values, bins=nbins, range=prange, normed=True,\
        xerrs=False, label="Control", color='k')
    he.set_label("ttW shape")
    ylim(1.e-4, .1)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    subplot(122)
    hist( a_sig.mctperp, weights=a_sig.weight, bins=nbins, range=prange, histtype="step", stacked=True,\
        normed=True, label="ttbar Shape", color='#440088')
    he = hist_errorbars( b_sig.mctperp.values, weights=b_sig.weight.values, bins=nbins, range=prange, normed=True,\
        xerrs=False, label="Control", color='k')
    he.set_label("ttW shape")
    ylim(0, 0.04)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")

    savefig("plots/top_ttw_shape{}.pdf".format(flavor))

def plot_ww_www_shape(flavor):
    """Compare the shapes of the wgamma and WZ backgrounds"""
    nbins=29
    prange=(10,300)

    a_sig = mc[smc['sig'+flavor]&(mc.mctype=="WWTo2L2Nu")]
    b_sig = mc[smc['sig'+flavor]&(mc.mctype=="WWW")]

    figure(figsize=(14,7))
    fig = subplot(121)
    fig.set_yscale('log', nonposy='clip')
    hist( a_sig.mctperp, weights=a_sig.weight, bins=nbins, range=prange, histtype="step", stacked=True,\
        normed=True, label="WW Shape", color='#440088')
    he = hist_errorbars( b_sig.mctperp.values, weights=b_sig.weight.values, bins=nbins, range=prange, normed=True,\
        xerrs=False, label="Control", color='k')
    he.set_label("WWW shape")
    ylim(1.e-4, .1)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    subplot(122)
    hist( a_sig.mctperp, weights=a_sig.weight, bins=nbins, range=prange, histtype="step", stacked=True,\
        normed=True, label="WW Shape", color='#440088')
    he = hist_errorbars( b_sig.mctperp.values, weights=b_sig.weight.values, bins=nbins, range=prange, normed=True,\
        xerrs=False, label="Control", color='k')
    he.set_label("WWW shape")
    ylim(0, 0.04)
    legend()
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")

    savefig("plots/ww_www_shape{}.pdf".format(flavor))
