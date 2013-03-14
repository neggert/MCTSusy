from config.data import *
import ROOT as R
import sys
sys.path.append("import")
import selection
from CMSPyLibs.plot import *
import numpy as np
from config.parameters import *

import CMSPyLibs.general_calc as general_calc

calc = general_calc.VarCalculator()

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('text.latex', preamble = '\usepackage{amsmath}')

from matplotlib.pyplot import *
from matplotlib.font_manager import *
switch_backend("pdf")

fontp = FontProperties(family="Helvetica", size=12)
fontpb = FontProperties(family="Helvetica", size=12, weight="book")

def recalc_MCT(event):
	p4jet1 = R.TLorentzVector(0., 0., 0., 0.)
	p4jet2 = R.TLorentzVector(0., 0., 0., 0.)
	p4l1 = R.TLorentzVector(0., 0., 0., 0.)
	p4l2 = R.TLorentzVector(0., 0., 0., 0.)
	p4l1.SetPtEtaPhiM(event['pt1'], event['eta1'], event['phi1'], 0.)
	p4l2.SetPtEtaPhiM(event['pt2'], event['eta2'], event['phi2'], 0.)

	# Add third lepton into MET
	# note that we don't care about eta or M for these vectors, as they're not used.
	p4met = R.TLorentzVector(0., 0., 0., 0.)
	p4met.SetPtEtaPhiM(event['metPt'], 0., event['metPhi'], 0.)
	p4l3 = R.TLorentzVector(0., 0., 0., 0.)
	p4l3.SetPtEtaPhiM(event['pt3'], 0., event['phi3'], 0.)

	p4met += p4l3

	calc.setP4s(p4jet1, p4jet2, p4l1, p4l2, p4met)
	return calc.mctPerp_210()

def recalc_MET(event):

	# Add third lepton into MET
	# note that we don't care about eta or M for these vectors, as they're not used.
	p4met = R.TLorentzVector(0., 0., 0., 0.)
	p4met.SetPtEtaPhiM(event['metPt'], 0., event['metPhi'], 0.)
	p4l3 = R.TLorentzVector(0., 0., 0., 0.)
	p4l3.SetPtEtaPhiM(event['pt3'], 0., event['phi3'], 0.)

	p4met += p4l3

	return p4met.Pt()

wz = mc[mc.ThirdLepton & ~np.isnan(mc.pt3) & (mc.metPt > 30)]

# do some shuffling to recover as many events as possible
def shuffle_leptons(event):
	# make an OF pair if possible
	if abs(event['pdg1']) == abs(event['pdg2']):
		if abs(event['pdg1']) != abs(event['pdg3']):
			event['pdg2'], event['pdg3'] = event['pdg3'], event['pdg2']
			event['pt2'], event['pt3'] = event['pt3'], event['pt2']
			event['eta2'], event['eta3'] = event['eta3'], event['eta2']
			event['phi2'], event['phi3'] = event['phi3'], event['phi2']

		elif abs(event['pdg2']) != abs(event['pdg3']):
			event['pdg1'], event['pdg3'] = event['pdg3'], event['pdg1']
			event['pt1'], event['pt3'] = event['pt3'], event['pt1']
			event['eta1'], event['eta3'] = event['eta3'], event['eta1']
			event['phi1'], event['phi3'] = event['phi3'], event['phi1']

	return event

wz = wz.apply(shuffle_leptons, axis=1)

# recalculate MCT

wz.mctperp = wz.apply(recalc_MCT, axis=1)
wz.mctperp *= 80.4/91.2
wz.metPt = wz.apply(recalc_MET, axis=1)
wz.ThirdLepton = False
swz = selection.get_samples(wz)
wz = wz[swz['sig_of']]

wz_data = data[data.ThirdLepton]
wz_data = wz_data.apply(shuffle_leptons, axis=1)
wz_data.mctperp = wz_data.apply(recalc_MCT, axis=1)
wz_data.mctperp *= 80.4/91.2
wz_data.metPt = wz_data.apply(recalc_MET, axis=1)
wz_data.ThirdLepton = False
swz_data = selection.get_samples(wz_data)
wz_data = wz_data[swz_data['sig_of']]


ww = mc[smc['sig_of'] & (mc.mc_cat=="WW")]

nbins=14
nrange=(10,290)

f = figure(figsize=(6,6))
f.set_facecolor('w')
fig = subplot2grid((4,1),(0,0), rowspan=3)
fig.set_yscale('log', nonposy='clip')
fig.set_ylim(0.001, 1000)
fig.set_ylabel("entries / 20 GeV", fontproperties=fontpb, color='k')

hist(ww.mctperp, weights=ww.weight, color=bkg_colors['WW'], bins=nbins, range=nrange, normed=True, histtype="step", label="WW MC")

he = hist_errorbars( wz.mctperp, weights=wz.weight, bins=nbins, range=nrange, normed=True,
    xerrs=False, color="k")

he.set_label(r"$3-1$ Lepton CR MC")
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

hist_ratio(wz.mctperp, ww.mctperp, ww.weight, wz.weight, bins=nbins, range=nrange, normed=True)
axhline(1, color="k")
fig2.set_ylim(0,2)
fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

figtext(0.12, 0.92, r"CMS Simulation $\sqrt{\text{s}}=8\;\text{TeV}$", color='k',
         fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

savefig("plots/3-1leptonCR.pdf")


#reweight
ww_mc, bins = np.histogram(ww.mctperp, weights=ww.weight, bins=nbins, range=nrange, normed=True)
wz_mc, bins = np.histogram(wz.mctperp, weights=wz.weight, bins=nbins, range=nrange, normed=True)
ratio = ww_mc/wz_mc

def get_weight(mctperp):
	i = np.argmax(np.where(bins<mctperp, bins, 0))
	return ratio[i]

reweighting = wz_data.mctperp.apply(get_weight)
wz_data.weight *= reweighting


f = figure(figsize=(6,6))
f.set_facecolor('w')
fig = subplot2grid((4,1),(0,0), rowspan=3)
fig.set_yscale('log', nonposy='clip')
fig.set_ylim(0.001, 1000)
fig.set_ylabel("entries / 20 GeV", fontproperties=fontpb, color='k')

hist(ww.mctperp, weights=ww.weight, color=bkg_colors['WW'], bins=nbins, range=nrange, normed=True, histtype="step", label="WW MC")

he = hist_errorbars( wz_data.mctperp, weights=wz_data.weight, bins=nbins, range=nrange, normed=True,
    xerrs=False, color="k")

he.set_label(r"$3-1$ Lepton Data")
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

hist_ratio(wz_data.mctperp, ww.mctperp, ww.weight, wz_data.weight, bins=nbins, range=nrange, normed=True)
axhline(1, color="k")
fig2.set_ylim(0,2)
fig2.set_ylabel("ratio", fontproperties=fontpb, color='k')

xlabel("$M_{\mathrm{CT}\perp}$ (GeV)", fontproperties=fontp, color='k')

figtext(0.12, 0.92, r"CMS Preliminary $\sqrt{\text{s}}=8\;\text{TeV},$\quad L$_{\text{int}}=19.5\;\text{fb}^{-1}$", color='k',
         fontproperties=FontProperties(family="Helvetica", size=12, weight="demi"))

savefig("plots/3-1leptonData.pdf")
