from pandas import *
import sys
sys.path.append("../")
import CMSPyLibs.plot
import selection
reload(CMSPyLibs.plot)
reload(selection)

import BackgroundFit
reload(BackgroundFit)

from data import *
from parameters import *

import validation
reload(validation)

import limits

from matplotlib.pyplot import *

channels =  ['_of', '_sf']

def make_validation_plots(flavor):
    validation.make_data_mc_plots(flavor)
    validation.plot_fake_closure(flavor)
    validation.plot_top_closure(flavor)
    validation.run_mc_closure(flavor)

def make_mc_plot(flavor):
    f = validation.plot_mc('sig'+flavor, 'mctperp', bins=30, plotrange=(0,300))
    f.set_yscale('log', nonposy='clip')
    f.set_ylim(0.01, 100000)
    xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
    savefig("plots/mc{}.pdf".format(flavor))

def bkg_predict(flavor):
    return BackgroundFit.do_bkg_fit(data, mc, 100., flavor)

def write_data_card(channels, outfile):
    f = open(outfile, 'w')
    f.write("lumi 9200 400\n")
    f.write("#<channel label> <observed> <total background> <total bkg uncertainty> <stat uncertainty> <wz> <wz uncertainty> <ttbar+fakes> <ttbar+fakes uncertainty> <zgamma> <zgamma uncertainty> <zz> <z uncertainty><rare> <rare uncertainty>\n")

    for i, ch in enumerate(channels):
        bkg_dict = bkg_predict(ch)
        sig = sum(data[sd['sig_mct_high'+ch]].weight)
        bkg, bkg_err = bkg_dict['high']['Total']
        output = " ".join(map(str, [i, sig, bkg, bkg_err, bkg_err, 0,0,0,0,0,0,0,0,0,0]))
        f.write(output)
        f.write("\n")

    f.close()

def write_sms_cards():
    limits.write_sms_dc(chi, "limits/8TeVc1c1.xsec", "/home/nic/cms/MCTSusy/2012HCP/limits/histo_chi.root", "limits/datacards/TChipmSlepSnu/TChipmSlepSnu", ["_of", "_sf"])
    limits.write_sms_dc(slep, "limits/8TeVeLeL.xsec", "/home/nic/cms/MCTSusy/2012HCP/limits/histo_slepslep.root", "limits/datacards/TSlepSlep/TSlepSlep", ["_sf",], xsec_scale=2.)

def write_data_cards():
    write_data_card(['_of', '_sf'], "limits/datacards/TChipmSlepSnu/data.txt")
    write_data_card(['_sf',], "limits/datacards/TSlepSlep/data.txt")



