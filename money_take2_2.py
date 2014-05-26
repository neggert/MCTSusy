import numpy as np
from config.parameters import *

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'weight':'bold'})
rc('text', usetex=True)
rc('text.latex', preamble = '\usepackage{amsmath}\usepackage{hepnicenames}\usepackage{hepunits}\usepackage{sansmath}\sansmath')
rc('xtick.major', size=7, width=1)
rc('xtick.minor', size=3, width=1)
rc('ytick.major', size=7, width=1)
rc('ytick.minor', size=3, width=1)
from matplotlib.ticker import MultipleLocator

import simplejson as json
import ROOT as R

from matplotlib.font_manager import *
import matplotlib.pyplot as plt
plt.switch_backend("Agg")


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

def extract_shape(ws, x, obs_name, shape_name, **kwargs):
    """Extract values from a RooFit PDF into an array"""
    obs = ws.obj(obs_name)
    shape = ws.obj(shape_name)
    out = np.zeros(x.shape)
    for i in xrange(len(x)):
        obs.setVal(x[i])
        out[i] = shape.getVal()
    return out

def extract_data_hist(ws, x, obs_name):
    data = ws.obj("obsData")
    out = np.zeros(x.shape)
    for i in xrange(data.numEntries()):
        data.get(i)
        out[i] = data.weight()
    return out

def extract_points_from_th1(th1, bins):
    y = []
    for i, b in enumerate(bins):
        y.append(th1.GetBinContent(i+1))

    return np.asarray(y)

sf_backgrounds = [
{"obs_name": "obs_x_sf", "shape_name": "wjets_sf_overallSyst_x_StatUncert_x_sf_wjets_syst_ShapeSys", "name":"wjets"},
{"obs_name": "obs_x_sf", "shape_name": "z_sf_overallSyst_x_StatUncert_x_sf_z_syst_ShapeSys", "name": "DY"},
{"obs_name": "obs_x_sf", "shape_name": "vv_sf_overallSyst_x_StatUncert", "name": "diboson"},
{"obs_name": "obs_x_sf", "shape_name": "of_sf_overallSyst_x_StatUncert", "name":"FS"}
]



def build_background_shape(ws, ch='sf', backgrounds=sf_backgrounds, log=True):

    with open("diboson_fracs2.json") as f:
        diboson_fracs = json.load(f)[ch]
    dibosons = ['ZZ', 'WZ']

    fontp = FontProperties(family="Helvetica", size=14)
    fontpb = FontProperties(family="Helvetica", size=12, weight="book")

    x = np.arange(10, 300, 10)
    bin_heights = [np.zeros(x.shape)]
    names = []
    for b in backgrounds:
        shape = extract_shape(ws, x=x, **b)
        if b['name'] == 'diboson':
            for process in dibosons:
                fracs = np.asarray(diboson_fracs[process])
                fracs[np.isnan(fracs) | np.isinf(fracs)] = 0.
                bin_heights.append(shape*fracs)
                names.append(process)
        else :
            bin_heights.append(shape)
            names.append(b['name'])

    fig = plt.figure(figsize=(6,6))
    fig.set_facecolor('w')
    ax = plt.subplot2grid((4,1),(0,0), rowspan=3)
    if log:
        ax.set_yscale('log', nonposy='clip')
        ax.set_ylim(0.01, 5000)
    else:
        ax.set_ylim(0., 1500)
    ax.set_ylabel("Entries / 10 GeV", fontproperties=fontp, color='k')


    top = bin_heights[0]
    for i in xrange(len(bin_heights)-1):
        bottom = top
        top = top+bin_heights[i+1]
        ax.fill(*get_step_fill_between(x, bottom, top), facecolor=bkg_colors[names[i]],
                label=bkg_labels[names[i]], closed=True)

    # plot SMS
    sms_file = R.TFile("signal_slep.root")
    sms1 = sms_file.Get("sms_template_{}_150_0".format(ch))
    sms2 = sms_file.Get("sms_template_{}_300_0".format(ch))

    sms1_points = extract_points_from_th1(sms1, x)
    sms2_points = extract_points_from_th1(sms2, x)

    sms_x = np.append(x, 300)
    sms1_points = np.append(sms1_points, sms1_points[-1])
    sms2_points = np.append(sms2_points, sms2_points[-1])


    ax.plot(sms_x, sms1_points, drawstyle="steps-post", linestyle="--", dashes=(3, 2), color="r", label=r'\noindent$m_{\PSlepton}=150\;\GeV$')
    ax.plot(sms_x, sms2_points, drawstyle="steps-post", linestyle="--", dashes=(3, 2), color="b", label=r'\noindent$m_{\PSlepton}=300\;\GeV$')

    data_height = extract_data_hist(ws, x, "obs_x_"+ch)
    data_x = x+5
    data_x_reduced = data_x[data_height>0]
    data_y = data_height[data_height > 0]
    data_err = np.sqrt(data_y)
    ax.errorbar(data_x_reduced, data_y, data_err, fmt='.', color='k', label="Data")


    # aesthetics
    ax.set_xlim(10, 300)
    handles, labels = ax.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()
    ax.legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=1, numpoints=1)
    ax.set_axisbelow(False)
    minorticks = MultipleLocator(10)
    ax.xaxis.set_minor_locator(minorticks)
    plt.setp(ax.get_xticklabels(), visible=False)

    # ratio plot
    ax2 = plt.subplot2grid((4,1),(3,0))
    xerrs = np.ones_like(data_x)*5

    for i in reversed(xrange(1, len(top))):
        if top[i] < 3:
            top[i-1] += top[i]
            top[i] = 0
            data_height[i-1] += data_height[i]
            data_height[i] = 0
            data_x[i-1] = data_x[i] - 5
            data_x[i] = 0
            xerrs[i-1] = xerrs[i] + 5
            xerrs[i] = 0

    vals = data_height/top
    errs = np.sqrt(data_height)/top
    ax2.errorbar(data_x[vals >0], vals[vals > 0], errs[vals>0], xerrs[vals > 0], color='k', fmt='.')
    plt.axhline(1, color="k")
    ax2.set_ylim(0.,2.0)
    ax2.set_xlim(10, 300)
    ax2.set_ylabel("Ratio", fontproperties=fontp, color='k')

    plt.xlabel("$M_{\mathrm{CT}\perp}$ [\GeV]", fontproperties=fontp, color='k')

    plt.figtext(0.5, 0.92, r"CMS \hspace{3em} $\sqrt{\text{s}}=8\;\TeV \hspace{3em} \text{L}=19.5\;\text{fb}^{-1}$", color='k', ha='center',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="bold"))
    chan_txt = r"$\Pepm\Pemp/\Pmupm\Pmump$"
    plt.figtext(0.42, 0.85, chan_txt, color='k', ha='center',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="bold"))

    if log:
        plt.savefig("plots/money{}2.pdf".format(ch))
    else:
        plt.savefig("plots/money{}2_linear.pdf".format(ch))



