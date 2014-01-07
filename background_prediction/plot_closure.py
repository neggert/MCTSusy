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


def get_step_fill_between(x, y1, y2):

    fill_x = np.zeros(2*len(x))
    fill_x[::2] = x
    fill_x[1:-1:2] = x[1:]
    fill_x[-1] = 300
    fill_y1 = np.zeros(fill_x.shape)
    fill_y1[::2] = y1
    fill_y1[1::2] = y1
    fill_y2 = np.zeros(fill_x.shape)
    fill_y2[::2] = y2
    fill_y2[1::2] = y2



    return np.append(fill_x, fill_x[::-1]), np.append(fill_y1, fill_y2[::-1])

def plot_top_closure(flavor):
    """Make closure plot for the top control sample"""
    truths = ['ttbar', 'tW', 'ttV', 'ttVV']
    truth = {}
    truth['ttbar'] = smc['sig_mct_low_'+flavor] & ((mc.mctype=='ttbar') | (mc.mctype=='TTG'))
    truth['tW'] = smc['sig_mct_low_'+flavor] & (mc.mctype=='tW')
    truth['ttV'] = smc['sig_mct_low_'+flavor] & ((mc.mctype=='ttw')|(mc.mctype=="TTZ"))
    truth ['ttVV'] = smc['sig_mct_low_'+flavor] & (mc.mctype=="TTWW")
    all_truth = smc['sig_mct_low_'+flavor] & (mc.mc_cat=='top')

    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.001, 10000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

    nbins = 19
    nrange = (10., 200.)
    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.001, 1000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

    x = []
    w = []
    for t in truths:
        vals = mc[truth[t]].mctperp.values
        vals[vals > nrange[1]] = nrange[1]-1e-6
        x.append(vals)
        w.append(mc[truth[t]].weight)

    c = [(.9,.6,0),
              (.5, .7, .9),
              (0,.6,.5),
              (.95, .9, .25)
             ]
    labels = [r't\=t', r'tW', r't\=tV', r't\=tVV']

    # hist( mc[toptruth].mctperp, weights=mc[toptruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
    #     normed=True, label="Top MC Truth", color=(.9,.6,0), linewidth=2)
    # hist( mc[twtruth].mctperp, weights=mc[twtruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
    #     normed=True, label="tW MC Truth", color=(0,.45,.7), linewidth=2)
    h1 = hist_errorbars(x, weights=w, color=c, bins=nbins, range=nrange, normed=True, stacked=True, histtype="step", label=labels, plotstyle="filled")

    values = mc[smc['top_mct_low_'+flavor]].mctperp.values
    values[values>nrange[1]] = nrange[1]-1e-6
    he = hist_errorbars( values, weights=mc[smc['top_mct_low_'+flavor]].weight.values, bins=nbins, range=nrange, normed=True,\
        xerrs=False, color="k")[2]

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

    v1 = mc[smc['top_ctrl_'+flavor]].mctperp.values
    v1[v1>nrange[1]] = nrange[1] -1.e-6
    v2 = mc[all_truth].mctperp.values
    v2[v2>nrange[1]] = nrange[1]-1.e-6

    hist_ratio(v1, v2, mc[all_truth].weight, mc[smc['top_ctrl_'+flavor]].weight, bins=nbins, range=nrange, normed=True)
    axhline(1, color="k")
    fig2.set_ylim(0,2)
    fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    figtext(0.12, 0.92, r"CMS Simulation $\sqrt{\text{s}}=8\;\text{TeV}$", color='k',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

    savefig("plots/closure_top_{}.pdf".format(flavor))

def print_top_closure():
    """Make closure plot for the top control sample"""

    names = ['top']

    truth = {}
    truth['top'] = mc[smc['sig_mct_low'] & (mc.mc_cat == 'top')]
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

    t = PrettyTable(["\mctp\ Cut",]+[str(c)+"\GeV" for c in cuts.tolist()])
    t.vertical_char = "&"
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
    wz_fs = (mc.mc_cat=="WZ") & ((abs(mc.parentParentPdg1) == 24) | (abs(mc.parentParentPdg1) == 24))
    truth['fs'] = mc[smc['sig_mct_low_sf'] & ((mc.mc_cat=='top') | (mc.mc_cat=='WW') | wz_fs | (mc.mc_cat=="VVV") | (mc.mc_cat=="HWW"))]
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

    t = PrettyTable(["\mctp\ Cut",]+[str(c)+"\GeV" for c in cuts.tolist()])
    t.vertical_char = "&"
    for i in xrange(len(names)):
        data_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], data[:,i])
        err_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], errs[:,i])
        combined = [datastr+" $\pm$ "+errstr for datastr, errstr in zip(data_strings, err_strings)]
        combined.insert(0, names[i])
        t.add_row(combined)
    data_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], control_data)
    err_strings = map(str.format, ["{:.5f}" for _ in xrange(len(cuts))], control_errs)
    combined = [datastr+" +- "+errstr for datastr, errstr in zip(data_strings, err_strings)]
    combined.insert(0, "CR")
    t.add_row(combined)

    print t

def plot_fake_closure():
    """Make closure plot for the top control sample"""
    truths = ['wjets', 'top']
    truth = {}
    truth['wjets'] = smc['sig_mct_low'] & (mc.mc_cat=="fake")
    truth['top'] = smc['sig_mct_low'] & (mc.mc_cat=='top') & (mc.gen_neutrinos <= 1)
    all_truth = truth['wjets'] | truth['top']


    nbins = 19
    nrange = (10., 200.)
    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.001, .1)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

    x = []
    w = []
    for t in truths:
        vals = mc[truth[t]].mctperp.values
        vals[vals > nrange[1]] = nrange[1]-1e-6
        x.append(vals)
        w.append(mc[truth[t]].weight)

    c = [bkg_colors['fake'], bkg_colors['top']
             ]
    labels = ['W+Jets', 'Non-dilepton Top']

    # hist( mc[toptruth].mctperp, weights=mc[toptruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
    #     normed=True, label="Top MC Truth", color=(.9,.6,0), linewidth=2)
    # hist( mc[twtruth].mctperp, weights=mc[twtruth].weight, bins=nbins, range=nrange, histtype="step", fill=False, rwidth=1.,\
    #     normed=True, label="tW MC Truth", color=(0,.45,.7), linewidth=2)
    h1 = hist_errorbars(x, weights=w, color=c, bins=nbins, range=nrange, normed=True, stacked=True, histtype="step", label=labels, plotstyle="filled")


    values = mc[smc['wjets_mct_low']].mctperp.values
    values[values>nrange[1]] = nrange[1]-1e-6
    he = hist_errorbars( values, weights=mc[smc['wjets_mct_low']].weight.values, bins=nbins, range=nrange, normed=True,\
        xerrs=False, color="k")[2]

    he.set_label("Control Region MC")
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

    v1 = mc[smc['wjets_ctrl']].mctperp.values
    v1[v1>nrange[1]] = nrange[1] -1.e-6
    v2 = mc[all_truth].mctperp.values
    v2[v2>nrange[1]] = nrange[1]-1.e-6

    hist_ratio(v1, v2, mc[all_truth].weight, mc[smc['wjets_ctrl']].weight, bins=nbins, range=nrange, normed=True)
    axhline(1, color="k")
    fig2.set_ylim(0,2)
    fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    figtext(0.12, 0.92, r"CMS Simulation $\sqrt{\text{s}}=8\;\text{TeV}$", color='k',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

    savefig("plots/closure_fake.pdf")


# def plot_fake_closure():
#     """Make closure plot for the fake lepton control sample"""
#     faketruth = smc['sig'] & (mc.mc_cat=='fake') | ((mc.mc_cat=="top") & (mc.gen_neutrinos==1))

#     f = figure(figsize=(6,6))
#     f.set_facecolor('w')
#     fig = subplot2grid((4,1),(0,0), rowspan=3)
#     fig.set_yscale('log', nonposy='clip')
#     fig.set_ylim(0.001, 10000)
#     fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

#     nbins = 9
#     nrange = (10., 100.)
#     hist( mc[faketruth].mctperp, weights=mc[faketruth].weight, bins=nbins, range=nrange, histtype="step", stacked=True,\
#         normed=True, label="W+Jets MC Truth", color=(.8,.4, 0), linewidth=2)
#     he = hist_errorbars( mc[smc['wjets_ctrl']].mctperp.values, weights=mc[smc['wjets_ctrl']].weight.values, bins=nbins, range=nrange, normed=True,\
#         xerrs=False, label="Control", color='k')
#     he.set_label("Control Region MC")
#     ylim(1.e-4, .5)

#     # move data to top of legend
#     handles, labels = fig.get_legend_handles_labels()
#     handles.insert(0,handles.pop())
#     labels.insert(0,labels.pop())

#     legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=1)
#     fig.set_axisbelow(False)

#     minorticks = MultipleLocator(10)
#     fig.xaxis.set_minor_locator(minorticks)

#     fig2 = subplot2grid((4,1),(3,0), sharex=fig)

#     hist_ratio(mc[smc['wjets_ctrl']].mctperp, mc[faketruth].mctperp, mc[faketruth].weight, mc[smc['wjets_ctrl']].weight, bins=nbins, range=nrange, normed=True)
#     axhline(1, color="k")
#     fig2.set_ylim(0,2)
#     fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

#     xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

#     figtext(0.12, 0.92, r"CMS Simulation $\sqrt{\text{s}}=8\;\text{TeV}$", color='k',
#              fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

#     savefig("plots/closure_fake.pdf")

def plot_fs_closure():
    """Make closure plot for the fake lepton control sample"""

    # the flavor symmetric part of WZ has one of the leptons from a W
    wz_fs = (mc.mc_cat=="WZ") & ((abs(mc.parentParentPdg1) == 24) | (abs(mc.parentParentPdg1) == 24))

    truths = ['top', 'WW', 'WZ', 'HWW', 'VVV']
    truth = {}
    truth['top'] = smc['sig_mct_low_sf'] & (mc.mc_cat=="top")
    truth['WW'] = smc['sig_mct_low_sf'] & (mc.mc_cat=="WW")
    truth['WZ'] = smc['sig_mct_low_sf'] & wz_fs
    truth['HWW'] = smc['sig_mct_low_sf'] & (mc.mc_cat=="HWW")
    truth['VVV'] = smc['sig_mct_low_sf'] & (mc.mc_cat=="VVV")

    fstruth = smc['sig_mct_low_sf'] & ((mc.mc_cat=='top') | (mc.mc_cat=='WW') | wz_fs | (mc.mc_cat=="VVV") | (mc.mc_cat=="HWW"))


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

    x = []
    w = []
    for t in truths:
        vals = mc[truth[t]].mctperp.values
        vals[vals > nrange[1]] = nrange[1]-1e-6
        x.append(vals)
        w.append(mc[truth[t]].weight)

    c = [(.9,.6,0),
          (.35, .7, .9),
          (0,.6,.5),
           (0, .45, .70),
           (.80,.40,0)]
    labels = [r'Top', r'WW', r'WZ', r'H$\rightarrow$WW', 'VVV']


    h1 = hist_errorbars(x, weights=w, color=c, bins=nbins, range=nrange, normed=True, stacked=True, histtype="step", label=labels, plotstyle="filled")
    values = mc[smc['sig_mct_low_of']].mctperp.values
    values[values>nrange[1]] = nrange[1]-1e-6
    he = hist_errorbars( values, weights=mc[smc['sig_mct_low_of']].weight.values, bins=nbins, range=nrange, normed=True,
        xerrs=False, label="Control", color='k')[2]
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

    v1 = mc[smc['sig_mct_low_of']].mctperp.values
    v1[v1>nrange[1]] = nrange[1] -1.e-6
    v2 = mc[fstruth].mctperp.values
    v2[v2>nrange[1]] = nrange[1]-1.e-6

    hist_ratio(v1, v2, mc[fstruth].weight,mc[smc['sig_mct_low_of']].weight,
               bins=nbins, range=nrange, normed=True)
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
