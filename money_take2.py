import numpy as np
from config.parameters import *

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('text.latex', preamble = '\usepackage{amsmath}')
from matplotlib.ticker import MultipleLocator

import simplejson as json

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
    is_sf = "sf" in obs_name
    data = ws.obj("obsData")
    if is_sf:
        data = data.reduce("channelCat==channelCat::sf")
    else:
        data = data.reduce("channelCat==channelCat::of")
    out = np.zeros(x.shape)
    for i in xrange(data.numEntries()):
        data.get(i)
        out[i] = data.weight()
    return out

sf_backgrounds = [
{"obs_name": "obs_x_sf", "shape_name": "wjets_sf_sf_overallSyst_x_StatUncert_x_sf_wjets_syst_sf_ShapeSys", "name":"wjets"},
{"obs_name": "obs_x_sf", "shape_name": "z_sf_sf_overallSyst_x_StatUncert_x_sf_z_syst_sf_ShapeSys", "name": "DY"},
{"obs_name": "obs_x_sf", "shape_name": "vv_sf_sf_overallSyst_x_StatUncert_x_sf_ww_syst_sf_ShapeSys", "name": "diboson"},
{"obs_name": "obs_x_sf", "shape_name": "top_sf_sf_overallSyst_x_StatUncert", "name":"top"}
]

of_backgrounds = [
{"obs_name": "obs_x_of", "shape_name": "wjets_of_of_overallSyst_x_StatUncert_x_of_wjets_syst_of_ShapeSys", "name":"wjets"},
{"obs_name": "obs_x_of", "shape_name": "z_of_of_overallSyst_x_StatUncert_x_of_z_syst_of_ShapeSys", "name": "DY"},
{"obs_name": "obs_x_of", "shape_name": "vv_of_of_overallSyst_x_StatUncert_x_of_ww_syst_of_ShapeSys", "name": "diboson"},
{"obs_name": "obs_x_of", "shape_name": "top_of_of_overallSyst_x_StatUncert", "name":"top"}
]


def build_background_shape(ws, ch='sf', backgrounds=sf_backgrounds, log=True):

    with open("diboson_fracs.json") as f:
        diboson_fracs = json.load(f)[ch]
    dibosons = ['Rare', 'ZZ', 'WZ', 'WW']

    fontp = FontProperties(family="Helvetica", size=12)
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
        ax.set_ylim(0.01, 2000)
    else:
        ax.set_ylim(0., 1500)
    ax.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')


    top = bin_heights[0]
    for i in xrange(len(bin_heights)-1):
        bottom = top
        top = top+bin_heights[i+1]
        ax.fill(*get_step_fill_between(x, bottom, top), facecolor=bkg_colors[names[i]],
                label=bkg_labels[names[i]], closed=True)

    data_height = extract_data_hist(ws, x, "obs_x"+ch)
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
    ax.legend(handles, labels, frameon=False, prop=fontpb, borderaxespad=1)
    ax.set_axisbelow(False)
    minorticks = MultipleLocator(10)
    ax.xaxis.set_minor_locator(minorticks)

    # ratio plot
    ax2 = plt.subplot2grid((4,1),(3,0))
    vals = data_height/top
    errs = np.sqrt(data_height)/top
    ax2.errorbar(data_x[vals >0], vals[vals > 0], errs[vals>0], color='k', fmt='.')
    plt.axhline(1, color="k")
    ax2.set_ylim(0.5,1.5)
    ax2.set_xlim(10, 300)
    ax2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    plt.xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    plt.figtext(0.12, 0.92, r"CMS Preliminary $\sqrt{\text{s}}=8\;\text{TeV},$\quad L$_{\text{int}}=19.5\;\text{fb}^{-1}$", color='k',
             fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

    if log:
        plt.savefig("plots/money{}.pdf".format(ch))
    else:
        plt.savefig("plots/money{}_linear.pdf".format(ch))



