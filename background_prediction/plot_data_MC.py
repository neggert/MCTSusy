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

def compare_data_mc(selection_name, variable, bins=20, plotrange=(0,100)):
    """
    Make a plot comparing the histogram of variable between data and simulation

    Arguments:
    selection_name -- name of the selection to be plotted. Should be one of those in the dict returned by get_samples in selection.pyplot
    variable -- variable to make a histogram of
    bins -- number of bins for the histogram (default 20)
    plotrange -- range to plot (default (0,100))
    """

    selected = mc[smc[selection_name]]
    data_selected = data[sd[selection_name]]

    # groups = selected.groupby('mc_cat')

    group_order = ['top', 'WV', 'ZZ', 'DY', 'fake']

    bkgtpl = []
    bkgwtpl = []
    bkgltpl = []

    for name in group_order:
        if selected[(selected.mc_cat==name) & (selected[variable] > plotrange[0])&(selected[variable] < plotrange[1])][variable].count() > 0:
            bkgtpl.append( selected[selected.mc_cat==name][variable])
            bkgwtpl.append( selected[selected.mc_cat==name].weight)
            bkgltpl.append(name)

    if len(bkgtpl) == 0:
        raise RuntimeError('No Events')

    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.001, 10000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

    h = hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl)
    print sum([sum(weights) for weights in bkgwtpl])

    he = hist_errorbars( data_selected[variable], xerrs=False, bins=bins, range=plotrange)
    he.set_label("Data")
    fig.set_axisbelow(False)

    # move data to top of legend
    handles, labels = fig.get_legend_handles_labels()
    handles.insert(0,handles.pop())
    labels.insert(0,labels.pop())

    legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=2)
    fig.set_axisbelow(False)

    minorticks = MultipleLocator(10)
    fig.xaxis.set_minor_locator(minorticks)

    fig2 = subplot2grid((4,1),(3,0), sharex=fig)
    hist_ratio(data_selected[variable], selected[variable], selected.weight, bins=bins, range=plotrange)

    axhline(1, color="k")
    fig2.set_ylim(0.5,1.5)
    fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    figtext(0.12, 0.92, r"CMS Preliminary $\sqrt{\text{s}}=8\;\text{TeV},$\quad L$_{\text{int}}=9.2\;\text{fb}^{-1}$", color='k',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

    return fig, fig2

def make_data_mc_plots( flavor):
    """Make all of the data-MC comparison plots for the given channel"""

    f,f2 = compare_data_mc('sig'+flavor, 'mctperp', 30, (0,300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 100000)
    f2.set_ylim(0.5, 1.5)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_sig{}.pdf".format(flavor))

    try:
        f,f2 = compare_data_mc('z_ctrl'+flavor, 'mctperp', 30, (0,300))
        f.set_yscale('log', nonposy='clip')
        f.set_ylim(0.01, 100000)
        f2.set_ylim(0.5, 1.5)
        xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
        savefig("plots/data_mc_z{}.pdf".format(flavor))

        f,f2 = compare_data_mc('z_ctrl'+flavor, 'mctperp', 19, (5,100))
        f.set_yscale('log', nonposy='clip')
        f.set_ylim(1, 1000)
        f2.set_ylabel("Data/MC")
        f2.set_ylim(0.5, 1.5)
        xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
        savefig("plots/data_mc_z_lowmct{}.pdf".format(flavor))
    except RuntimeError:
        pass

    try :
        f,f2 = compare_data_mc('wjets_ctrl'+flavor, 'mctperp', 10, (0,100))
        f.set_yscale('log', nonposy='clip')
        f.set_ylim(0.01, 100000)
        f2.set_ylim(0.5, 1.5)
        xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
        savefig("plots/data_mc_fake{}.pdf".format(flavor))
    except RuntimeError:
        pass

    f,f2 = compare_data_mc('top_ctrl'+flavor, 'mctperp', 20, (0,200))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 100000)
    f2.set_ylim(0.5, 1.5)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_top{}.pdf".format(flavor))

    f,f2 = compare_data_mc('3lep_ctrl'+flavor, 'mctperp', 20, (0,200))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 100000)
    f2.set_ylim(0.5, 1.5)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_3l{}.pdf".format(flavor))

    f,f2 = compare_data_mc('1tag_ctrl'+flavor, 'mctperp', 20, (0,200))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 100000)
    f2.set_ylim(0.5, 1.5)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_1tag{}.pdf".format(flavor))

    f,f2 = compare_data_mc('2tag_ctrl'+flavor, 'mctperp', 20, (0,200))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 100000)
    f2.set_ylim(0.5, 1.5)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_2tag{}.pdf".format(flavor))

    f,f2 = compare_data_mc('moretag_ctrl'+flavor, 'mctperp', 20, (0,200))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 100000)
    f2.set_ylim(0.5, 1.5)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_moretag{}.pdf".format(flavor))

def plot_mc(selection_name, variable, bins=20, plotrange=(0,100)):
    selected = mc[smc[selection_name]]

    # groups = selected.groupby('mc_cat')

    group_order = ['top', 'WV', 'ZZ', 'DY', 'fake']

    bkgtpl = []
    bkgwtpl = []
    bkgltpl = []

    for name in group_order:
        if selected[(selected.mc_cat==name) & (selected.mctperp > plotrange[0])&(selected.mctperp < plotrange[1])].mctperp.count() > 0:
            bkgtpl.append( selected[selected.mc_cat==name][variable])
            bkgwtpl.append( selected[selected.mc_cat==name].weight)
            bkgltpl.append(name)

    if len(bkgtpl) == 0:
        raise RuntimeError('No Events')

    figure()
    fig = subplot(111)
    h = hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl)
    fig.set_axisbelow(False)

    legend()

    return fig

# if __name__ == '__main__':
#     make_data_mc_plots('')


