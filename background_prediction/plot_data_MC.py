import sys
sys.path.append("import/")
import CMSPyLibs.plot
reload(CMSPyLibs.plot)
# from CMSPyLibs.plot import *
from config.data import *
from config.parameters import bkg_colors, bkg_labels, lumi
from prettytable import PrettyTable


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'weight':'bold'})
rc('text', usetex=True)
rc('text.latex', preamble = '\usepackage{amsmath}\usepackage{hepnicenames}\usepackage{hepunits}\usepackage{sansmath}\sansmath')
rc('xtick.major', size=7, width=1)
rc('xtick.minor', size=3, width=1)
rc('ytick.major', size=7, width=1)
rc('ytick.minor', size=3, width=1)

from matplotlib.pyplot import *
from matplotlib.font_manager import *
switch_backend("pdf")

fontp = FontProperties(family="Helvetica", size=12)
fontpb = FontProperties(family="Helvetica", size=10, weight="book")

def load_xsec(filename):
    """
    load cross-sections for a file. Returns a dictionary describing cross-sections on their uncertainties.
    The dictionary is indexed by the particle mass, and each element is a tuple containing (xsec, uncertainty)
    """
    f = open(filename)
    xsec_dict = {}
    for line in f:
        mass, xsec, err = line.split()
        xsec_dict[float(mass)] = (float(xsec), float(err))
    f.close()
    return xsec_dict

def compare_data_mc(selection_name, variable, bins=20, plotrange=(0,100), cumulative=False):
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

    selected.mc_cat[selected.mc_cat.isin(["VVV", "HWW"])] = "Rare"

    # groups = selected.groupby('mc_cat')

    group_order = ['top', 'WW', 'WZ', 'ZZ', 'Rare', 'DY', 'fake']
    group_order.reverse()

    bkgtpl = []
    bkgwtpl = []
    bkgltpl = []
    bkgctpl = []

    for name in group_order:
        if selected[(selected.mc_cat==name) & (selected[variable] > plotrange[0])][variable].count() > 0:
            vars = selected[selected.mc_cat==name][variable].values
            vars[vars>plotrange[1]] = plotrange[1]-1e-6

            bkgtpl.append( selected[selected.mc_cat==name][variable])
            bkgwtpl.append( selected[selected.mc_cat==name].weight)
            bkgltpl.append(bkg_labels[name])
            bkgctpl.append(bkg_colors[name])

    if len(bkgtpl) == 0:
        raise RuntimeError('No Events')

    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    # fig.set_yscale('log', nonposy='clip')
    # fig.set_ylim(0.001, 10000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

    h = CMSPyLibs.plot.hist_errorbars(bkgtpl, plotstyle="filled", weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl, color=bkgctpl, linewidth=0.5)
    print sum([sum(weights) for weights in bkgwtpl])

    datavars =  data_selected[variable].values
    datavars[datavars>plotrange[1]] = plotrange[1]-1.e-6
    he = CMSPyLibs.plot.hist_errorbars( data_selected[variable], xerrs=False, bins=bins, range=plotrange)
    he[-1].set_label("Data")
    fig.set_axisbelow(False)

    # move data to top of legend
    handles, labels = fig.get_legend_handles_labels()
    handles.insert(0,handles.pop())
    labels.insert(0,labels.pop())

    l = legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=1)
    fig.set_axisbelow(False)

    minorticks = MultipleLocator(10)
    fig.xaxis.set_minor_locator(minorticks)

    fig2 = subplot2grid((4,1),(3,0), sharex=fig)
    CMSPyLibs.plot.hist_ratio(data_selected[variable], selected[variable], selected.weight, bins=bins, range=plotrange)

    axhline(1, color="k")
    fig2.set_ylim(0.5,1.5)
    fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    figtext(0.5, 0.92, r"CMS Unpublished \hspace{3em} $\sqrt{\text{s}}=8\;\TeV \hspace{3em} \text{L}=19.5\;\text{fb}^{-1}$", color='k', ha='center',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="bold"))

    return fig, fig2

def make_data_mc_plots():
    """Make all of the data-MC comparison plots for the given channel"""

    channels = ['of', 'sf']

    for flavor in channels:
        f,f2 = compare_data_mc('sig_'+flavor, 'mctperp', 29, (10, 300))
        f.set_yscale('log', nonposy='clip')
        f.set_ylim(0.01, 100000)
        f2.set_ylim(0, 2)
        xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
        savefig("plots/data_mc_sig_{}.pdf".format(flavor))

        f,f2 = compare_data_mc('top_ctrl_'+flavor, 'mctperp', 29, (10, 300))
        f.set_yscale('log', nonposy='clip')
        f.set_ylim(0.01, 100000)
        f2.set_ylim(0, 2)
        xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
        savefig("plots/data_mc_top_{}.pdf".format(flavor))

    f,f2 = compare_data_mc('z_ctrl_sf', 'mctperp', 29, (10, 300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 10000)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_z.pdf")

    f,f2 = compare_data_mc('z_ctrl_0met', 'mctperp', 29, (10, 300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 1000000)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_z_0met.pdf")

    f,f2 = compare_data_mc('z_ctrl_0met', 'metPt', 30, (0,300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(.1, 1e7)
    f2.set_ylim(0, 2)
    xlabel("MET (GeV)")
    savefig("plots/data_mc_z_metdist.pdf")

    f,f2 = compare_data_mc('z_ctrl_0met', 'mll', 30, (76,106))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(1e4, 1e7)
    f2.set_ylim(0, 2)
    xlabel("MET (GeV)")
    savefig("plots/data_mc_z_mll.pdf")

    f,f2 = compare_data_mc('z_ctrl_30met', 'mctperp', 29, (10, 300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 1000000)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_z_30met.pdf")

    f,f2 = compare_data_mc('z_ctrl_sf', 'mctperp', 9, (10,100))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(1, 10000)
    f2.set_ylabel("Data/MC")
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_z_lowmct.pdf")

    f,f2 = compare_data_mc('wjets_ctrl_sf', 'mctperp', 9, (10,100))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 100000)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_fake.pdf")

    f,f2 = compare_data_mc('wz_ctrl', 'mctperp', 29, (10, 300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 100)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_3l.pdf")

    f = plot_mc("sig", 'mctperp', 29, (10, 300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 10000)
    f.set_xlim(10, 300)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/mc_only.pdf")

    f = plot_mc("sig_of", 'mctperp', 29, (10, 300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 5000)
    f.set_xlim(10, 300)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/mc_only_of.pdf")

    f = plot_mc("sig_sf", 'mctperp', 29, (10, 300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 5000)
    f.set_xlim(10, 300)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/mc_only_sf.pdf")

    for flavor in channels:
        f,f2 = compare_data_mc('sig_'+flavor, 'mctperp', 29, (10, 300))
        # f.set_ylim(0.01, 5000)
        f2.set_ylim(0, 2)
        xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
        savefig("plots/data_mc_sig_{}_linear.pdf".format(flavor))

        f,f2 = compare_data_mc('top_ctrl_'+flavor, 'mctperp', 29, (10, 300))
        # f.set_ylim(0.01, 100000)
        f2.set_ylim(0, 2)
        xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
        savefig("plots/data_mc_top_{}_linear.pdf".format(flavor))

    f,f2 = compare_data_mc('z_ctrl_sf', 'mctperp', 29, (10, 300))
    # f.set_ylim(0.01, 10000)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_z_linear.pdf")

    f,f2 = compare_data_mc('z_ctrl_0met', 'mctperp', 29, (10, 300))
    # f.set_ylim(0.01, 1000000)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_z_0met_linear.pdf")

    f,f2 = compare_data_mc('z_ctrl_0met', 'metPt', 30, (0,300))
    # f.set_ylim(.1, 1e7)
    f2.set_ylim(0, 2)
    xlabel("MET (GeV)")
    savefig("plots/data_mc_z_metdist_linear.pdf")

    f,f2 = compare_data_mc('z_ctrl_30met', 'mctperp', 29, (10, 300))
    # f.set_ylim(0.01, 1000000)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_z_30met_linear.pdf")

    f,f2 = compare_data_mc('z_ctrl_sf', 'mctperp', 9, (10,100))
    # f.set_ylim(1, 10000)
    f2.set_ylabel("Data/MC")
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_z_lowmct_linear.pdf")

    f,f2 = compare_data_mc('wjets_ctrl_sf', 'mctperp', 9, (10,100))
    # f.set_ylim(0.01, 100000)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_fake_linear.pdf")

    f,f2 = compare_data_mc('wz_ctrl', 'mctperp', 29, (10, 300))
    # f.set_ylim(0.01, 100)
    f2.set_ylim(0, 2)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/data_mc_3l_linear.pdf")

    f = plot_mc("sig", 'mctperp', 29, (10, 300))
    # f.set_ylim(0.01, 10000)
    f.set_xlim(10, 300)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/mc_only_linear.pdf")

    f = plot_mc("sig_of", 'mctperp', 29, (10, 300))
    # f.set_ylim(0.01, 5000)
    f.set_xlim(10, 300)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/mc_only_of_linear.pdf")

    f = plot_mc("sig_sf", 'mctperp', 29, (10, 300))
    # f.set_ylim(0.01, 5000)
    f.set_xlim(10, 300)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/mc_only_sf_linear.pdf")

def plot_mc(selection_name, variable, bins=20, plotrange=(0,100)):
    selected = mc[smc[selection_name]]

    selected.mc_cat[selected.mc_cat.isin(["VVV", "HWW"])] = "Rare"

    # groups = selected.groupby('mc_cat')

    group_order = ['top', 'WW', 'WZ', 'ZZ', 'Rare', 'DY', 'fake']
    group_order.reverse()

    bkgtpl = []
    bkgwtpl = []
    bkgltpl = []
    bkgctpl = []

    for name in group_order:
        if selected[(selected.mc_cat==name) & (selected[variable] > plotrange[0])&(selected[variable] < plotrange[1])][variable].count() > 0:
            bkgtpl.append( selected[selected.mc_cat==name][variable])
            bkgwtpl.append( selected[selected.mc_cat==name].weight)
            bkgltpl.append(bkg_labels[name])
            bkgctpl.append(bkg_colors[name])

    if len(bkgtpl) == 0:
        raise RuntimeError('No Events')

    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot(111)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.01, 50000)
    fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')

    h = hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl, color=bkgctpl, linewidth=0.5)
    print sum([sum(weights) for weights in bkgwtpl])

    # plot some example signals on top


    hist_filename = "limits/TChipmSlepSnu_Nevents.root"
    xsec_filename = "limits/8TeVc1c1.xsec"
    if not os.path.exists(hist_filename):
        raise IOError(hist_filename+" does not exist.")
    nevents_file = R.TFile(hist_filename)
    nevents_hist = nevents_file.Get("ScanValidator/NEvents_histo")
    xsec_dict = load_xsec(xsec_filename)

    sms1 = chi[sel_chi[selection_name] & (chi.mass1==200) & (chi.mass2==25)]
    sms1.weight *= xsec_dict[200][0] / nevents_hist.GetBinContent(nevents_hist.FindBin(200, 25)) * lumi
    sms3 = chi[sel_chi[selection_name] & (chi.mass1==400) & (chi.mass2==100)]
    sms3.weight *= xsec_dict[400][0] / nevents_hist.GetBinContent(nevents_hist.FindBin(400, 100)) * lumi


    x = [sms1.mctperp, sms3.mctperp]
    w = [sms1.weight, sms3.weight]
    c = ['b','r']
    labels = [r'\noindent$m_{\chi^\pm}=200$ GeV, $m_{\chi^0}=25$ GeV', r'\noindent$m_{\chi^\pm}=400$ GeV, $m_{\chi^0}=100$ GeV']

    print h[0][-1]
    print hist(x, histtype="step", weights=w, color=c, label=labels, bins=h[1], bottom=h[0][-1], ls='dashed', zorder=10, lw=2)
    # m, bins = np.histogram(sms1.mctperp, weights=sms1.mctperp, bins=h[1])
    # m += h[0][-1]
    # # draw by hand
    # x = np.zeros(2 * len(bins) - 1, np.float)
    # x[0:2*len(bins)-1:2], x[1:2*len(bins)-1:2] = bins, bins[:-1]
    # y = np.zeros(2 * len(bins) - 1, np.float)
    # y[1:2*len(bins)-1:2], y[2:2*len(bins):2] = m, m
    # fill(x, y, closed=False, edgecolor="b", fill=False)

    # print x, y

    handles, labels = fig.get_legend_handles_labels()

    fig.set_xlim(*plotrange)

    legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=1, ncol=2)
    fig.set_axisbelow(False)

    minorticks = MultipleLocator(10)
    fig.xaxis.set_minor_locator(minorticks)

    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    figtext(0.12, 0.92, r"CMS Simulation $\sqrt{\text{s}}=8\;\text{TeV},$\quad L$_{\text{int}}=19.5\;\text{fb}^{-1}$", color='k',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

    return fig

def print_ZZ_datamc():
    """Make closure plot for the top control sample"""
    flavor = 'sf'
    names = ['top', 'tW']

    from config.data import *

    truth = {}
    mc_cr = mc[smc['z_mct_low_'+flavor]]
    data_cr = data[sd['z_mct_low_'+flavor]]

    cuts = np.arange(100, 320, 20)

    data = np.zeros(len(cuts))
    errs = np.zeros(data.shape)
    control_data = np.zeros(len(cuts))
    control_errs = np.zeros(control_data.shape)

    for i, cut in enumerate(cuts):
        data[i] = mc_cr[mc_cr.mctperp > cut].weight.sum()
        errs[i] = np.sqrt(sum(mc_cr[mc_cr.mctperp > cut].weight**2))
        control_data[i] = data_cr[data_cr.mctperp > cut].weight.sum()
        control_errs[i] = np.sqrt(sum(data_cr[data_cr.mctperp > cut].weight**2))

    t = PrettyTable(["",]+[str(c)+"\GeV" for c in cuts.tolist()])
    t.vertical_char = "&"
    data_strings = map(str.format, ["{:.1f}" for _ in xrange(len(cuts))], data)
    err_strings = map(str.format, ["{:.1f}" for _ in xrange(len(cuts))], errs)
    combined = [datastr+" $\pm$ "+errstr for datastr, errstr in zip(data_strings, err_strings)]
    combined.insert(0, "Monte Carlo")
    t.add_row(combined)
    data_strings = map(str.format, ["{:d}" for _ in xrange(len(cuts))], map(int, control_data))
    err_strings = map(str.format, ["{:f}" for _ in xrange(len(cuts))], control_errs)
    combined = data_strings
    combined.insert(0, "Data")
    t.add_row(combined)

    print t

def print_WZ_datamc():
    """Make closure plot for the top control sample"""
    flavor = 'sf'
    names = ['top', 'tW']

    from config.data import *

    truth = {}
    mc_cr = mc[smc['wz_ctrl']]
    data_cr = data[sd['wz_ctrl']]

    cuts = np.arange(20, 300, 20)

    data = np.zeros(len(cuts))
    errs = np.zeros(data.shape)
    control_data = np.zeros(len(cuts))
    control_errs = np.zeros(control_data.shape)

    for i, cut in enumerate(cuts):
        data[i] = mc_cr[mc_cr.mctperp > cut].weight.sum()
        errs[i] = np.sqrt(sum(mc_cr[mc_cr.mctperp > cut].weight**2))
        control_data[i] = data_cr[data_cr.mctperp > cut].weight.sum()
        control_errs[i] = np.sqrt(sum(data_cr[data_cr.mctperp > cut].weight**2))

    t = PrettyTable(["\mctp\ Cut",]+map(str, cuts.tolist()))
    t.vertical_char = "&"
    data_strings = map(str.format, ["{:.1f}" for _ in xrange(len(cuts))], data)
    err_strings = map(str.format, ["{:.1f}" for _ in xrange(len(cuts))], errs)
    combined = [datastr+" $\pm$ "+errstr for datastr, errstr in zip(data_strings, err_strings)]
    combined.insert(0, "Monte Carlo")
    t.add_row(combined)
    data_strings = map(str.format, ["{:d}" for _ in xrange(len(cuts))], map(int, control_data))
    err_strings = map(str.format, ["{:f}" for _ in xrange(len(cuts))], control_errs)
    combined = data_strings
    combined.insert(0, "Data")
    t.add_row(combined)

    print t



if __name__ == '__main__':
    # make_data_mc_plots()
    f = plot_mc('sig_mct_low_sf', 'mctperp', 29, (10,300))
    savefig("plots/mc_only_sf.pdf")
    f = plot_mc('sig_mct_low_of', 'mctperp', 29, (10,300))
    savefig("plots/mc_only_of.pdf")



