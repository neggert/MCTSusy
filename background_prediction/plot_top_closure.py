import sys
sys.path.append("../../")
sys.path.append("import/")
from CMSPyLibs.plot import *
from config.data import *


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('text.latex', preamble = '\usepackage{amsmath}')

from matplotlib.pyplot import *
from matplotlib.font_manager import *
switch_backend("pdf")

fontp = FontProperties(family="Helvetica", size=12)
fontpb = FontProperties(family="Helvetica", size=12, weight="book")


def plot_top_closure(flavor):
    """Make closure plot for the top control sample"""
    toptruth = smc['sig_'+flavor] & (mc.mc_cat == 'top') &(mc.mctype!='tW')
    twtruth = smc['sig_'+flavor] & (mc.mctype=='tW')

    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.001, 10000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

    nbins = 29
    nrange = (10., 300.)
    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.001, 1000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')
    hist( mc[toptruth].mctperp, weights=mc[toptruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
        normed=True, label="Top MC Truth", color=(.9,.6,0), linewidth=2)
    hist( mc[twtruth].mctperp, weights=mc[twtruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
        normed=True, label="tW MC Truth", color=(0,.45,.7), linewidth=2)
    he = hist_errorbars( mc[smc['top_ctrl_'+flavor]].mctperp.values, weights=mc[smc['top_ctrl_'+flavor]].weight.values, bins=nbins, range=nrange, normed=True,\
        xerrs=False, color="k")
    he.set_label("$\geq$1 b-tags MC")
    ylim(1.e-7, .1)
    fig.set_axisbelow(False)

    # move data to top of legend
    handles, labels = fig.get_legend_handles_labels()
    handles.insert(0,handles.pop())
    labels.insert(0,labels.pop())

    legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=1)
    fig.set_axisbelow(False)

    minorticks = MultipleLocator(10)
    fig.xaxis.set_minor_locator(minorticks)

    fig2 = subplot2grid((4,1),(3,0), sharex=fig)

    hist_ratio(mc[smc['top_ctrl_'+flavor]].mctperp, mc[toptruth | twtruth].mctperp, mc[toptruth | twtruth].weight, bins=nbins, range=nrange, normed=True)
    axhline(1, color="k")
    fig2.set_ylim(0,2)
    fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    figtext(0.12, 0.92, r"CMS Simulation $\sqrt{\text{s}}=8\;\text{TeV}$", color='k',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

    savefig("plots/closure_top_{}.pdf".format(flavor))

def plot_fake_closure():
    """Make closure plot for the fake lepton control sample"""
    faketruth = smc['sig'] & (mc.mc_cat=='wjets')

    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.001, 10000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

    nbins = 9
    nrange = (10., 100.)
    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.001, 1000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')
    hist( mc[faketruth].mctperp, weights=mc[faketruth].weight, bins=nbins, range=nrange, histtype="step", stacked=True,\
        normed=True, label="W+Jets MC Truth", color=(.8,.4, 0), linewidth=2)
    he = hist_errorbars( mc[smc['wjets_ctrl']].mctperp.values, weights=mc[smc['wjets_ctrl']].weight.values, bins=nbins, range=nrange, normed=True,\
        xerrs=False, label="Control", color='k')
    he.set_label("Control Region MC")
    ylim(1.e-4, .5)

    # move data to top of legend
    handles, labels = fig.get_legend_handles_labels()
    handles.insert(0,handles.pop())
    labels.insert(0,labels.pop())

    legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=1)
    fig.set_axisbelow(False)

    minorticks = MultipleLocator(10)
    fig.xaxis.set_minor_locator(minorticks)

    fig2 = subplot2grid((4,1),(3,0), sharex=fig)

    hist_ratio(mc[smc['wjets_ctrl']].mctperp, mc[faketruth].mctperp, mc[faketruth].weight, bins=nbins, range=nrange, normed=True)
    axhline(1, color="k")
    fig2.set_ylim(0,2)
    fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    figtext(0.12, 0.92, r"CMS Simulation $\sqrt{\text{s}}=8\;\text{TeV}$", color='k',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

    savefig("plots/closure_fake.pdf")


if __name__ == '__main__':
    plot_top_closure('of')
    plot_top_closure('sf')
    plot_fake_closure()
