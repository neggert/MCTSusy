#!/usr/bin/env python

import ROOT as R
import sys
import re

f = sys.argv[1]
rfile = R.TFile(f)


ws = rfile.Get("combined")

data = ws.data("obsData")

sbmodel = ws.obj("ModelConfig")

constr = sbmodel.GetNuisanceParameters()
R.RooStats.RemoveConstantParameters(constr)

# run the initial fit
sbmodel.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Minos(sbmodel.GetParametersOfInterest()))

# make fit result dictionary

results = {'sf':{'top': (sbmodel.GetNuisanceParameters().find("n_sf_top").getVal(), 0),
		         'vv': (sbmodel.GetNuisanceParameters().find("n_sf_vv").getVal(), 0),
		         'DY': (sbmodel.GetNuisanceParameters().find("n_sf_z").getVal(), 0),
		         'fake': (sbmodel.GetNuisanceParameters().find("n_sf_wjets").getVal(), 0)
		        },
		   'of':{'top': (sbmodel.GetNuisanceParameters().find("n_of_top").getVal(), 0),
				 'vv': (sbmodel.GetNuisanceParameters().find("n_of_vv").getVal(), 0),
				 'DY': (sbmodel.GetNuisanceParameters().find("n_of_z").getVal(), 0),
				 'fake': (sbmodel.GetNuisanceParameters().find("n_of_wjets").getVal(), 0)
				}
}

poi_hat = sbmodel.GetParametersOfInterest().first().getVal()
poi_hat_err = sbmodel.GetParametersOfInterest().first().getError()

print poi_hat, poi_hat_err

from config.data import *
from prep_hists import load_xsec

m1, m2 = map(int, re.search("_(\d+)_(\d+)_", f).groups())

if not (re.match(".*slep.*", f) is None):
	sms = slep
	sel = sel_slep 
	xsec_file = "limits/8TeVeLeL.xsec"
	xsec_multiplier = 2.
	hist_filename = "limits/TSlepSlep_Nevents.root"


elif not (re.match(".*chi.*", f) is None):
	sms = chi
	sel = sel_chi
	xsec_file = "limits/8TeVc1c1.xsec"
	xsec_multiplier = 1.
	hist_filename = "limits/TChipmSlepSnu_Nevents.root"

nevents_file = R.TFile(hist_filename)
nevents_hist = nevents_file.Get("ScanValidator/NEvents_histo")
events_per_point = nevents_hist.GetBinContent(nevents_hist.FindBin(m1, m2))

xsec_dict = load_xsec(xsec_file)
xsec, _ = xsec_dict[float(m1)]
xsec *= xsec_multiplier

sms.weight *= xsec*lumi/events_per_point

from matplotlib.pyplot import *
sys.path.append("import/")
sys.path.append("../import/")
from config.data import *
from config.parameters import bkg_colors
from CMSPyLibs.plot import *

bins=29
plotrange = (10,300)
figs = []
for ch in ['of','sf']:
    bkgtpl = []
    bkgwtpl = []
    bkgltpl = []
    bkgctpl = []

    bkgtpl.append( sms[sel['sig_mct_low_'+ch] & (sms.mass1==m1) & (sms.mass2 == m2)].mctperp )

    bkgwtpl.append( poi_hat*sms[sel['sig_mct_low_'+ch] & (sms.mass1==m1) & (sms.mass2 == m2)].weight )
    bkgltpl.append("SMS")
    bkgctpl.append('r')

    bkgtpl.append( data[sd['top_ctrl_'+ch]].mctperp )
    data_norm = data[sd['top_mct_low_'+ch]].mctperp.count()
    est_events = float(results[ch]['top'][0])
    sf = est_events/data_norm
    bkgwtpl.append( sf*np.ones(data[sd['top_ctrl_'+ch]].mctperp.count()) )
    bkgltpl.append("Top")
    bkgctpl.append(bkg_colors['top'])

    bkgtpl.append( mc[smc['sig_'+ch]&(mc.mc_cat=="WW")].mctperp )
    allvv_norm = sum(mc[smc['sig_mct_low_'+ch]&((mc.mc_cat=='WW') | (mc.mc_cat=='ZZ') | (mc.mc_cat=="WV"))].weight)
    est_events = float(results[ch]['vv'][0])
    sf = est_events/allvv_norm
    bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mc_cat=="WW")].weight)
    bkgltpl.append("WW")
    bkgctpl.append(bkg_colors['WW'])

    bkgtpl.append( mc[smc['sig_'+ch]&(mc.mc_cat=="WZ")].mctperp )
    bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mc_cat=="WZ")].weight)
    bkgltpl.append("WZ")
    bkgctpl.append(bkg_colors['WZ'])

    bkgtpl.append( mc[smc['sig_'+ch]&(mc.mc_cat=="ZZ")].mctperp )
    bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mc_cat=="ZZ")].weight)
    bkgltpl.append("ZZ")
    bkgctpl.append(bkg_colors['ZZ'])

    bkgtpl.append( mc[smc['sig_'+ch]&(mc.mc_cat=="VVV")].mctperp )
    bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mc_cat=="VVV")].weight)
    bkgltpl.append("VVV")
    bkgctpl.append(bkg_colors['VVV'])

    bkgtpl.append( mc[smc['sig_'+ch]&(mc.mc_cat=="HWW")].mctperp )
    bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mc_cat=="HWW")].weight)
    bkgltpl.append("HWW")
    bkgctpl.append(bkg_colors['HWW'])

    bkgtpl.append( mc[smc['sig_'+ch]&(mc.mc_cat=="DY")].mctperp )
    data_norm = sum(mc[smc['sig_mct_low_'+ch]&(mc.mc_cat=="DY")].weight)
    est_events = float(results[ch]['DY'][0])
    sf = est_events/data_norm
    bkgwtpl.append( sf*mc[smc['sig_'+ch]&(mc.mc_cat=="DY")].weight )
    bkgltpl.append("Z/$\gamma^*$")
    bkgctpl.append(bkg_colors['Z'])

    bkgtpl.append( data[sd['wjets_ctrl_'+ch]].mctperp )
    data_norm = data[sd['wjets_mct_low_'+ch]].mctperp.count()
    est_events = float(results[ch]['fake'][0])
    sf = est_events/data_norm
    bkgwtpl.append( sf*np.ones(data[sd['wjets_ctrl_'+ch]].mctperp.count()) )
    bkgltpl.append("Non-prompt")
    bkgctpl.append(bkg_colors['wjets'])

    f = figure(figsize=(6,6))
    f.set_facecolor('w')
    fig = subplot2grid((4,1),(0,0), rowspan=3)
    fig.set_yscale('log', nonposy='clip')
    fig.set_ylim(0.01, 2000)
    # fig.set_ylabel("entries / 10 GeV", fontproperties=fontpb, color='k')
    h = hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl,
             zorder=1, linewidth=0.5, color=bkgctpl)
    he = hist_errorbars( data[sd['sig_'+ch]].mctperp, xerrs=False, bins=bins, range=plotrange)
    he.set_label("Data")

    # move data to top of legend
    handles, labels = fig.get_legend_handles_labels()
    handles.insert(0,handles.pop())
    labels.insert(0,labels.pop())

    legend(handles, labels, frameon=False, borderaxespad=1)
    fig.set_axisbelow(False)

    minorticks = MultipleLocator(10)
    fig.xaxis.set_minor_locator(minorticks)

    all_bkg = Series()
    all_bkg_w = Series()
    for mct, w in zip(bkgtpl, bkgwtpl):
        all_bkg = all_bkg.append(mct)
        all_bkg_w = all_bkg_w.append(Series(w))


    fig2 = subplot2grid((4,1),(3,0))
    hist_ratio(data[sd['sig_'+ch]].mctperp, all_bkg, all_bkg_w, bins=bins, range=plotrange)
    axhline(1, color="k")
    fig2.set_ylim(0.0,2.0)
    # fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

    # xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

    # figtext(0.12, 0.92, r"CMS Preliminary $\sqrt{\text{s}}=8\;\text{TeV},$\quad L$_{\text{int}}=19.5\;\text{fb}^{-1}$", color='k',
             # fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))


show()
raw_input("...")
