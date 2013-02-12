import sys
sys.path.append("../../")
sys.path.append("import/")
from CMSPyLibs.plot import *
from config.data import *
from config.parameters import bkg_colors
from prettytable import PrettyTable
import numpy as np


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

def print_top_closure():
    """Make closure plot for the top control sample"""

    names = ['top', 'tW']

    truth = {}
    truth['top'] = mc[smc['sig_mct_low'] & (mc.mc_cat == 'top') &(mc.mctype!='tW')]
    truth['tW'] = mc[smc['sig_mct_low'] & (mc.mctype=='tW')]
    control = mc[smc['top_mct_low']]

    cuts = np.arange(20, 220, 20)

    data = np.zeros((len(cuts), len(names)))
    errs = np.zeros(data.shape)
    control_data = np.zeros(len(cuts))
    control_errs = np.zeros(control_data.shape)

    for i, cut in enumerate(cuts):
        for j, name in enumerate(names):
            data[i, j] = truth[name][truth[name].mctperp > cut].weight.sum()/truth[name].weight.sum()
            errs[i,j ] = np.sqrt(sum(truth[name][truth[name].mctperp > cut].weight**2))/truth[name].weight.sum()
        control_data[i] = control[control.mctperp > cut].weight.sum()/control.weight.sum()
        control_errs[i] = np.sqrt(sum(control[control.mctperp > cut].weight**2))/control.weight.sum()

    t = PrettyTable(["",]+map(str, cuts.tolist()))
    for i in xrange(len(names)):
        data_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], data[:,i])
        err_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], errs[:,i])
        combined = [datastr+" +- "+errstr for datastr, errstr in zip(data_strings, err_strings)]
        combined.insert(0, names[i])
        t.add_row(combined)
    data_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], control_data)
    err_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], control_errs)
    combined = [datastr+" +- "+errstr for datastr, errstr in zip(data_strings, err_strings)]
    combined.insert(0, "CR")
    t.add_row(combined)

    print t

def print_fs_closure():
    """Make closure plot for the top control sample"""

    names = ['fs']

    truth = {}
    truth['fs'] = mc[smc['sig_mct_low_sf'] & ((mc.mc_cat=='top') | (mc.mc_cat=='WW') | (mc.mc_cat=="WZ"))]
    control = mc[smc['sig_mct_low_of']]

    cuts = np.arange(20, 260, 20)

    data = np.zeros((len(cuts), len(names)))
    errs = np.zeros(data.shape)
    control_data = np.zeros(len(cuts))
    control_errs = np.zeros(control_data.shape)

    for i, cut in enumerate(cuts):
        for j, name in enumerate(names):
            data[i, j] = truth[name][truth[name].mctperp > cut].weight.sum()/truth[name].weight.sum()
            errs[i,j ] = np.sqrt(sum(truth[name][truth[name].mctperp > cut].weight**2))/truth[name].weight.sum()
        control_data[i] = control[control.mctperp > cut].weight.sum()/control.weight.sum()
        control_errs[i] = np.sqrt(sum(control[control.mctperp > cut].weight**2))/control.weight.sum()

    t = PrettyTable(["",]+map(str, cuts.tolist()))
    for i in xrange(len(names)):
        data_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], data[:,i])
        err_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], errs[:,i])
        combined = [datastr+" +- "+errstr for datastr, errstr in zip(data_strings, err_strings)]
        combined.insert(0, names[i])
        t.add_row(combined)
    data_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], control_data)
    err_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], control_errs)
    combined = [datastr+" +- "+errstr for datastr, errstr in zip(data_strings, err_strings)]
    combined.insert(0, "CR")
    t.add_row(combined)

    print t

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

def plot_fs_closure():
    """Make closure plot for the fake lepton control sample"""
    fstruth = smc['sig_sf'] & ((mc.mc_cat=='top') | (mc.mc_cat=='WW') | (mc.mc_cat=="WZ"))

    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.001, 10000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

    nbins = 29
    nrange = (10., 300.)

    hist( mc[fstruth].mctperp, weights=mc[fstruth].weight, bins=nbins, range=nrange, histtype="step", stacked=True,\
        normed=True, label="Flavor Symmetric MC Truth", color=bkg_colors['WW'], linewidth=2)
    he = hist_errorbars( mc[smc['sig_of']].mctperp.values, weights=mc[smc['sig_of']].weight.values, bins=nbins, range=nrange, normed=True,\
        xerrs=False, label="Control", color='k')
    he.set_label("Control Region MC")
    ylim(1.e-7, .1)

    # move data to top of legend
    handles, labels = fig.get_legend_handles_labels()
    handles.insert(0,handles.pop())
    labels.insert(0,labels.pop())

    legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=1)
    fig.set_axisbelow(False)

    minorticks = MultipleLocator(10)
    fig.xaxis.set_minor_locator(minorticks)

    fig2 = subplot2grid((4,1),(3,0), sharex=fig)

    hist_ratio(mc[smc['sig_of']].mctperp, mc[fstruth].mctperp, mc[fstruth].weight, bins=nbins, range=nrange, normed=True)
    axhline(1, color="k")
    fig2.set_ylim(0,2)
    fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    figtext(0.12, 0.92, r"CMS Simulation $\sqrt{\text{s}}=8\;\text{TeV}$", color='k',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

    savefig("plots/closure_fs.pdf")


if __name__ == '__main__':
    plot_top_closure('of')
    plot_top_closure('sf')
    plot_fake_closure()
    plot_fs_closure()
