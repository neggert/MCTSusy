# -*- coding: utf-8 -*-
# <nbformat>3</nbformat>

# <markdowncell>

# Background $M_{CT}$ Shapes
# ===================
# The idea here is to find the shapes of the $M_{CT}$ distributions for different backgrounds using control samples.
# 
# Control regions:
# 
# *    top: Require both leptons to pass isolation cuts, but invert the b-tag veto, requireing $\geq$ 1 b-tag. This control region models the ttbar and tW backgrounds
# *    wjets:  B-tag veto, but require 1 of the leptons to be in the isolation control region, while the other lepton passes the tight isolation cut. This control region models the W+jets background
# *    WV: The shape of the WW/WZ background is determined from WW Monte Carlo
# *    ZZ: The shape is taken from Monte Carlo

# <markdowncell>

# Load libraries and the data

# <codecell>

from pandas import *
import scipy.stats
import sys
sys.path.append("../")
from CMSPyLibs.plot import *
from selection import *
from rootutils import *

# <codecell>

import ROOT as r
r.gSystem.SetDynamicPath(r.gSystem.Getenv("LD_LIBRARY_PATH"));
r.gSystem.Load("libRooFitCore")
r.gSystem.Load("libRooFit")
r.gSystem.Load("libRooStats")

# <codecell>

store = HDFStore("Data/background.hdf5")
data1 = store['ttbar']
data2 = store['tW']
data3 = store['WW'][:store['WW'].weight.count()/2]
data3.weight *= 2
data4 = store['wjets']
data5 = store['DY'][:store['DY'].weight.count()/2]
data5.weight *= 2
data6 = store['WZ']
data7 = store['ZZ'][:store['ZZ'].weight.count()/2]
data7.weight *= 2
data = data1.append(data2, ignore_index=True)
data = data.append(data3, ignore_index=True)
data = data.append(data4, ignore_index=True)
data = data.append(data5, ignore_index=True)
data = data.append(data6, ignore_index=True)
data = data.append(data7, ignore_index=True)
mcww = store['WW'][store['WW'].weight.count()/2:]
mcww.weight *= 2
mcz = store['DY'][store['DY'].weight.count()/2:]
mcz.weight *= 2
mczz = store['ZZ'][store['ZZ'].weight.count()/2:]
mczz.weight *= 2
mctw = store['tW']
mcwz = store['WZ']

# <markdowncell>

# Run selections on the data. This creates all the control and signal samples. Results are stored in sel. More detailed explanations of the different selections are in the docstring for get_samples.
# 
# We can now get whatever information we need from the control and signal samples

# <codecell>

cut = 100.
sel = get_samples( data, cut)
selww = get_samples( mcww, cut)
selzz = get_samples( mczz, cut)
selz = get_samples( mcz, cut)
seltw = get_samples( mctw, cut)
selwz = get_samples( mcwz, cut)
print sel.keys()

# <markdowncell>

# Get the true number of events from each mc type in the low $M_{\text{CT}}$ control region with b-jet veto and isolation cuts applied

# <codecell>

ntop_low_mct_true =  sum(data[sel['sig_mct_low'] & (data.mctype=="ttbar")].weight) + sum(data[sel['sig_mct_low'] & (data.mctype=="tW")].weight)
nwjets_low_mct_true =  sum(data[sel['sig_mct_low'] & (data.mctype=="wjets")].weight)
nww_low_mct_true =  sum(data[sel['sig_mct_low'] & ((data.mctype=="WW") | (data.mctype=="WZ"))].weight)
nzz_low_mct_true =  sum(data[sel['sig_mct_low'] & (data.mctype=="ZZ")].weight)
nz_low_mct_true =  sum(data[sel['sig_mct_low'] & (data.mctype=="DY")].weight)

# <markdowncell>

# Background Normalization
# ========================
# 
# Given a shape for a background, all we need to make a prediction for the number of events in the high signal region is some sort of normalization. We get the normalization by fitting for the number of events from each background in the low $M_{CT}$ control region, but which pass all other cuts.
# 
# Control region:
# 
# *    low $M_{CT}$: Apply b-tag veto and lepton isolation cuts. $5 < M_{CT} < 100$ GeV
# 
# The extended PDF for this region is given by adding up the shapes from the control region. Here $S_x(M_{CT})$ indicates the shape as a function of $M_{CT}$ as obtained in control region x
# 
# $$ P(M_{CT}) = N_{\text{top}}S_{\text{top}}(M_{CT}) + N_{\text{wjets}}S_{\text{wjets}}(M_{CT}) + N_{\text{WW}}S_{\text{WW}}(M_{CT}) $$
# 
# We use an extended likelihood fit to find $N_{\text{top}}$, $N_{\text{wjets}}$, and $N_{\text{WW}}$.

# <markdowncell>

# Load the data from the low $M_{CT}$ control region into a RooDataSet

# <codecell>

mct = r.RooRealVar("mct", "mct", 5., 100.)
w = r.RooRealVar("w", "w", 0., 10.)
ds = create_roodataset( data[sel['sig_mct_low']].mctperp, data[sel['sig_mct_low']].weight, mct, w, title="data")

# <markdowncell>

# Build the background shapes using events in the control regions

# <codecell>

# first make a TH1
topshapehistTH1 = create_TH1( data[sel['top_mct_low']].mctperp, data[sel['top_mct_low']].weight, "top")
wjetsshapehistTH1 = create_TH1( data[sel['wjets_mct_low']].mctperp, data[sel['wjets_mct_low']].weight, "wjets")
wwshapehistTH1 = create_TH1( mcww[selww['sig_mct_low']].mctperp, mcww[selww['sig_mct_low']].weight, "ww")
zzshapehistTH1 = create_TH1( mczz[selzz['sig_mct_low']].mctperp, mczz[selzz['sig_mct_low']].weight, "zz")
zshapehistTH1 = create_TH1( mcz[selz['sig_mct_low']].mctperp, mcz[selz['sig_mct_low']].weight, "Z")
# Then convert to a RooDataHist
topshapehist = r.RooDataHist("topHist", "topHist", r.RooArgList(mct), topshapehistTH1)
wjetsshapehist = r.RooDataHist("wjetsHist", "wjetsHist", r.RooArgList(mct), wjetsshapehistTH1)
wwshapehist = r.RooDataHist("wwHst", "wwHist", r.RooArgList(mct), wwshapehistTH1)
zzshapehist = r.RooDataHist("zzHst", "zzHist", r.RooArgList(mct), zzshapehistTH1)
zshapehist = r.RooDataHist("zHist", "zHist", r.RooArgList(mct), zshapehistTH1)
# and create a PDF from the RooDataHist
topshape = r.RooHistPdf( "top", "top", r.RooArgSet(mct), topshapehist)
wjetsshape = r.RooHistPdf( "wjets", "wjets", r.RooArgSet(mct), wjetsshapehist)
wwshape = r.RooHistPdf( "ww", "ww", r.RooArgSet(mct), wwshapehist)
zzshape = r.RooHistPdf( "zz", "zz", r.RooArgSet(mct), zzshapehist)
zshape = r.RooHistPdf("z", "z", r.RooArgSet(mct), zshapehist)

# <codecell>

# build the extended PDFs
ntop = r.RooRealVar("ntop", "ntop", 514, -10000., 10000.)
topext = r.RooExtendPdf( "topext", "topext", topshape, ntop)
nwjets = r.RooRealVar("nwjets", "nwjets", 300, -10000., 10000.)
wjetsext = r.RooExtendPdf( "wjetsext", "wjetsext", wjetsshape, nwjets)
nww = r.RooRealVar("nww", "nww", 600, -10000., 10000.)
wwext = r.RooExtendPdf( "wwext", "ww", wwshape, nww)
#nzz = r.RooFormulaVar("nzz", "nww*"+str(sum(mczz[selzz['sig_mct_low']].weight)/(sum(mcww[selww['sig_mct_low']].weight)+sum(mcwz[selwz['sig_mct_low']].weight))), r.RooArgList(nww))
nzz = r.RooRealVar("nzz", "nzz", 600, -10000., 10000.)
zzext = r.RooExtendPdf( "zzext", "zz", zzshape, nzz)
nz = r.RooRealVar("nz", "nz", 300, -10000, 10000)
zext = r.RooExtendPdf("zext", "z", zshape, nz)

# <markdowncell>

# Add up the shapes to create a composite model

# <codecell>

model = r.RooAddPdf("model", "model", r.RooArgList(topext, wjetsext, wwext, zzext, zext))

# <markdowncell>

# Using b-tagging efficiency to constrain $N_\mathrm{top}$
# ========================================================
# Because the top and WW shapes are very similar in the low $M_{CT}$ region, the fit tends to have very large uncertainties on these numbers. We can improve this by providing an indpendent estimate of $N_{\text{top}}$ and putting it in the fit as a Gaussian constraint.
# 
# First we estimate the b-tagging efficiency. Assume that any event that has at least 1 b-tag actually has 2 b-tags. Let $N_x$ be the number of events with exactly $x$ b-tags and $\epsilon$ be the b-tagging efficiency. Then
# 
# $$ N_2 = N_{\text{bb}}\epsilon^2 $$
# $$ N_1 = 2N_{\text{bb}}\epsilon(1-\epsilon) $$
# $$ N_0 = N_{\text{bb}}(1-\epsilon)^2 $$
# 
# Then we get
# 
# $$ \epsilon = \frac{2N_2}{N_1+2N_2} $$
# 
# with uncertainty
# 
# $$ \sigma_\epsilon = \epsilon\sqrt{\frac{1}{2N_2}-\frac{1}{N_1+2N_2}} $$

# <codecell>

n_1tag_tt = sum(data[sel['1tag_mct_low']].weight) - sum(mctw[seltw['1tag_mct_low']].weight)
n_1tag = sum(data[sel['1tag_mct_low']].weight)
n_2tag = sum(data[sel['2tag_mct_low']].weight)
eff = 2.*n_2tag/(n_1tag_tt+2*n_2tag)
seff = eff*sqrt(1./2/n_2tag - 1./(n_1tag_tt+2*n_2tag))
print eff, "+-", seff

# <markdowncell>

# We can estimate the number of bb+X events with 0 b-tags, which is the number of tt/tW events in our low $M_{CT}$ control region
# 
# $$N_0 = \frac{N_1}{2}\frac{1-\epsilon}{\epsilon} $$
# 
# with uncertainty
# 
# $$\sigma_{N_0} = N_0\sqrt{\frac{2}{N_1}+4\frac{\sigma_\epsilon^2}{\epsilon^2}} $$

# <codecell>

ntop_pred = n_1tag/2.*(1-eff)/eff
ntop_pred_err = ntop_pred*sqrt(2./n_1tag+4*seff**2/eff**2)
print "Actual:", ntop_low_mct_true
print "Predicted:", ntop_pred, "+-", ntop_pred_err

# <markdowncell>

# Create the Gaussian constraint

# <codecell>

# make a Gaussian constraint on nww
topMean = r.RooRealVar("topMean","topMean", ntop_pred)
topSigma = r.RooRealVar("topSigma","topSigma", ntop_pred_err)
topConstraint = r.RooGaussian("topConstraint","topContraint",ntop,topMean,topSigma)

# <codecell>

ntop.setVal(topMean.getVal())
ntop.setError(topSigma.getVal())
ntop.setConstant(False)

# <markdowncell>

# ##Constrain ratio of WW+WZ to ZZ using MC

# <codecell>

nwv_to_nzz_ratio = sum(mczz[selzz['sig_mct_low']].weight)/(sum(mcww[selww['sig_mct_low']].weight)+sum(mcwz[selwz['sig_mct_low']].weight))
zzMean = r.RooFormulaVar("zzMean", "nww*"+str(nwv_to_nzz_ratio), r.RooArgList(nww))
zzSigma = r.RooFormulaVar("zzSigma", "sqrt(nww)*"+str(nwv_to_nzz_ratio), r.RooArgList(nww))
zzConstraint = r.RooGaussian("zzConstraint", "zzConstraint", nzz, zzMean, zzSigma)
nzz.setVal(zzMean.getVal())
nzz.setError(zzSigma.getVal())
nzz.setConstant(False)

# <markdowncell>

# Set up the other parameters as constraints so they'll be included in the fit result

# <codecell>

ntop_ctrl_high_val = sum(data[sel['top_mct_high']].weight)
ntop_ctrl_high_err_val = scipy.stats.poisson.std(ntop_ctrl_high_val)
ntop_ctrl_high_var = r.RooRealVar("ntop_ctrl_high", "ntop_ctrl_high", ntop_ctrl_high_val, 0, 10000)
ntop_ctrl_high_obs = r.RooRealVar("ntop_ctrl_high_obs", "ntop_ctrl_high_obs", ntop_ctrl_high_val)
ntop_ctrl_high_err = r.RooRealVar("ntop_ctrl_high_err", "ntop_ctrl_high_err", ntop_ctrl_high_err_val)
ntop_ctrl_high_constraint = r.RooGaussian("ntop_ctrl_high_constraint", "ntop_ctrl_high_constraint", ntop_ctrl_high_var,
                                           ntop_ctrl_high_obs, ntop_ctrl_high_err)
ntop_ctrl_low_val = sum(data[sel['top_mct_low']].weight)
ntop_ctrl_low_err_val = scipy.stats.poisson.std(ntop_ctrl_low_val)
ntop_ctrl_low_var = r.RooRealVar("ntop_ctrl_low", "ntop_ctrl_low", ntop_ctrl_low_val, 0, 10000)
ntop_ctrl_low_obs = r.RooRealVar("ntop_ctrl_low_obs", "ntop_ctrl_low_obs", ntop_ctrl_low_val)
ntop_ctrl_low_err = r.RooRealVar("ntop_ctrl_low_err", "ntop_ctrl_low_err", ntop_ctrl_low_err_val)
ntop_ctrl_low_constraint = r.RooGaussian("ntop_ctrl_low_constraint", "ntop_ctrl_low_constraint", ntop_ctrl_low_var,
                                           ntop_ctrl_low_obs, ntop_ctrl_low_err)
nz_ctrl_high_val = sum(mcz[selz['sig_mct_high']].weight)
nz_ctrl_high_err_val = scipy.stats.poisson.std(nz_ctrl_high_val) if nz_ctrl_high_val > 0 else 1
nz_ctrl_high_var = r.RooRealVar("nz_ctrl_high", "nz_ctrl_high", nz_ctrl_high_val, 0, 10000)
nz_ctrl_high_obs = r.RooRealVar("nz_ctrl_high_obs", "nz_ctrl_high_obs", nz_ctrl_high_val)
nz_ctrl_high_err = r.RooRealVar("nz_ctrl_high_err", "nz_ctrl_high_err", nz_ctrl_high_err_val)
nz_ctrl_high_constraint = r.RooGaussian("nz_ctrl_high_constraint", "nz_ctrl_high_constraint", nz_ctrl_high_var,
                                           nz_ctrl_high_obs, nz_ctrl_high_err)
nz_ctrl_low_val = sum(mcz[selz['sig_mct_low']].weight)
nz_ctrl_low_err_val = scipy.stats.poisson.std(nz_ctrl_low_val)
nz_ctrl_low_var = r.RooRealVar("nz_ctrl_low", "nz_ctrl_low", nz_ctrl_low_val, 0, 10000)
nz_ctrl_low_obs = r.RooRealVar("nz_ctrl_low_obs", "nz_ctrl_low_obs", nz_ctrl_low_val)
nz_ctrl_low_err = r.RooRealVar("nz_ctrl_low_err", "nz_ctrl_low_err", nz_ctrl_low_err_val)
nz_ctrl_low_constraint = r.RooGaussian("nz_ctrl_low_constraint", "nz_ctrl_low_constraint", nz_ctrl_low_var,
                                           nz_ctrl_low_obs, nz_ctrl_low_err)
nwjets_ctrl_high_val = sum(data[sel['wjets_mct_high']].weight)
nwjets_ctrl_high_err_val = scipy.stats.poisson.std(nwjets_ctrl_high_val)
nwjets_ctrl_high_var = r.RooRealVar("nwjets_ctrl_high", "nwjets_ctrl_high", nwjets_ctrl_high_val, 0, 10000)
nwjets_ctrl_high_obs = r.RooRealVar("nwjets_ctrl_high_obs", "nwjets_ctrl_high_obs", nwjets_ctrl_high_val)
nwjets_ctrl_high_err = r.RooRealVar("nwjets_ctrl_high_err", "nwjets_ctrl_high_err", nwjets_ctrl_high_err_val)
nwjets_ctrl_high_constraint = r.RooGaussian("nwjets_ctrl_high_constraint", "nwjets_ctrl_high_constraint", nwjets_ctrl_high_var,
                                           nwjets_ctrl_high_obs, nwjets_ctrl_high_err)
nwjets_ctrl_low_val = sum(data[sel['wjets_mct_low']].weight)
nwjets_ctrl_low_err_val = scipy.stats.poisson.std(nwjets_ctrl_low_val)
nwjets_ctrl_low_var = r.RooRealVar("nwjets_ctrl_low", "nwjets_ctrl_low", nwjets_ctrl_low_val, 0, 10000)
nwjets_ctrl_low_obs = r.RooRealVar("nwjets_ctrl_low_obs", "nwjets_ctrl_low_obs", nwjets_ctrl_low_val)
nwjets_ctrl_low_err = r.RooRealVar("nwjets_ctrl_low_err", "nwjets_ctrl_low_err", nwjets_ctrl_low_err_val)
nwjets_ctrl_low_constraint = r.RooGaussian("nwjets_ctrl_low_constraint", "nwjets_ctrl_low_constraint", nwjets_ctrl_low_var,
                                           nwjets_ctrl_low_obs, nwjets_ctrl_low_err)
nww_ctrl_high_val = sum(mcww[selww['sig_mct_high']].weight)
nww_ctrl_high_err_val = scipy.stats.poisson.std(nww_ctrl_high_val)
nww_ctrl_high_var = r.RooRealVar("nww_ctrl_high", "nww_ctrl_high", nww_ctrl_high_val, 0, 10000)
nww_ctrl_high_obs = r.RooRealVar("nww_ctrl_high_obs", "nww_ctrl_high_obs", nww_ctrl_high_val)
nww_ctrl_high_err = r.RooRealVar("nww_ctrl_high_err", "nww_ctrl_high_err", nww_ctrl_high_err_val)
nww_ctrl_high_constraint = r.RooGaussian("nww_ctrl_high_constraint", "nww_ctrl_high_constraint", nww_ctrl_high_var,
                                           nww_ctrl_high_obs, nww_ctrl_high_err)
nww_ctrl_low_val = sum(mcww[selww['sig_mct_low']].weight)
nww_ctrl_low_err_val = scipy.stats.poisson.std(nww_ctrl_low_val)
nww_ctrl_low_var = r.RooRealVar("nww_ctrl_low", "nww_ctrl_low", nww_ctrl_low_val, 0, 10000)
nww_ctrl_low_obs = r.RooRealVar("nww_ctrl_low_obs", "nww_ctrl_low_obs", nww_ctrl_low_val)
nww_ctrl_low_err = r.RooRealVar("nww_ctrl_low_err", "nww_ctrl_low_err", nww_ctrl_low_err_val)
nww_ctrl_low_constraint = r.RooGaussian("nww_ctrl_low_constraint", "nww_ctrl_low_constraint", nww_ctrl_low_var,
                                           nww_ctrl_low_obs, nww_ctrl_low_err)
nzz_ctrl_high_val = sum(mczz[selzz['sig_mct_high']].weight)
nzz_ctrl_high_err_val = scipy.stats.poisson.std(nzz_ctrl_high_val)
nzz_ctrl_high_var = r.RooRealVar("nzz_ctrl_high", "nzz_ctrl_high", nzz_ctrl_high_val, 0, 10000)
nzz_ctrl_high_obs = r.RooRealVar("nzz_ctrl_high_obs", "nzz_ctrl_high_obs", nzz_ctrl_high_val)
nzz_ctrl_high_err = r.RooRealVar("nzz_ctrl_high_err", "nzz_ctrl_high_err", nzz_ctrl_high_err_val)
nzz_ctrl_high_constraint = r.RooGaussian("nzz_ctrl_high_constraint", "nzz_ctrl_high_constraint", nzz_ctrl_high_var,
                                           nzz_ctrl_high_obs, nzz_ctrl_high_err)
nzz_ctrl_low_val = sum(mczz[selzz['sig_mct_low']].weight)
nzz_ctrl_low_err_val = scipy.stats.poisson.std(nzz_ctrl_low_val)
nzz_ctrl_low_var = r.RooRealVar("nzz_ctrl_low", "nzz_ctrl_low", nzz_ctrl_low_val, 0, 10000)
nzz_ctrl_low_obs = r.RooRealVar("nzz_ctrl_low_obs", "nzz_ctrl_low_obs", nzz_ctrl_low_val)
nzz_ctrl_low_err = r.RooRealVar("nzz_ctrl_low_err", "nzz_ctrl_low_err", nzz_ctrl_low_err_val)
nzz_ctrl_low_constraint = r.RooGaussian("nzz_ctrl_low_constraint", "nzz_ctrl_low_constraint", nzz_ctrl_low_var,
                                           nzz_ctrl_low_obs, nzz_ctrl_low_err)
print nz_ctrl_high_val, nz_ctrl_high_err_val
print nwjets_ctrl_high_val, nwjets_ctrl_low_val

# <codecell>

all_pdfs = r.RooArgList(ntop_ctrl_low_constraint, ntop_ctrl_high_constraint, nz_ctrl_low_constraint, nz_ctrl_high_constraint,
                               nwjets_ctrl_low_constraint, nwjets_ctrl_high_constraint, nww_ctrl_low_constraint, nww_ctrl_high_constraint)
all_pdfs.add(r.RooArgList(nzz_ctrl_low_constraint, nzz_ctrl_high_constraint, model))
constraint_vars = r.RooArgSet(ntop_ctrl_low_var, ntop_ctrl_high_var, nz_ctrl_low_var, nz_ctrl_high_var,
                              nwjets_ctrl_low_var, nwjets_ctrl_high_var, nww_ctrl_low_var, nww_ctrl_high_var)
constraint_vars.add(r.RooArgSet(nzz_ctrl_low_var, nzz_ctrl_high_var))

# <codecell>

model_c = r.RooProdPdf("model_c", "model_c", all_pdfs)

# <markdowncell>

# Do the fit

# <codecell>

results = model_c.fitTo( ds, r.RooFit.ExternalConstraints(r.RooArgSet(topConstraint, zzConstraint)),
    r.RooFit.Constrain(constraint_vars), r.RooFit.Save())

# <markdowncell>

# Set up all of the calculated variables for the actual background yields. Doing it this way lets RooFit automagically calculate the errors for us, including the full covariance matrix.

# <codecell>

ntop_sig_pred = r.RooFormulaVar("ntop_sig_pred", "ntop_sig_pred", "ntop*ntop_ctrl_high/ntop_ctrl_low",
                                r.RooArgList(ntop, ntop_ctrl_high_var, ntop_ctrl_low_var))
nz_sig_pred = r.RooFormulaVar("nz_sig_pred", "nz_sig_pred", "nz*nz_ctrl_high/nz_ctrl_low",
                                r.RooArgList(nz, nz_ctrl_high_var, nz_ctrl_low_var))
nww_sig_pred = r.RooFormulaVar("nww_sig_pred", "nww_sig_pred", "nww*nww_ctrl_high/nww_ctrl_low",
                                r.RooArgList(nww, nww_ctrl_high_var, nww_ctrl_low_var))
nzz_sig_pred = r.RooFormulaVar("nzz_sig_pred", "nzz_sig_pred", "nzz*nzz_ctrl_high/nzz_ctrl_low",
                                r.RooArgList(nzz, nzz_ctrl_high_var, nzz_ctrl_low_var))
nwjets_sig_pred = r.RooFormulaVar("nwjets_sig_pred", "nwjets_sig_pred", "nwjets*nwjets_ctrl_high/nwjets_ctrl_low",
                                r.RooArgList(nwjets, nwjets_ctrl_high_var, nwjets_ctrl_low_var))
ntot_sig_pred = r.RooFormulaVar("ntot_sig_pred", "ntot_sig_pred", "ntop_sig_pred+nz_sig_pred+nww_sig_pred+nwjets_sig_pred+nzz_sig_pred",
                                r.RooArgList(ntop_sig_pred, nz_sig_pred, nww_sig_pred, nwjets_sig_pred, nzz_sig_pred))

# <markdowncell>

# Plot

# <codecell>

r.gROOT.SetStyle("Plain")

f = mct.frame()
ds.plotOn(f, r.RooFit.Binning(19, 5., 100.))
model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected))
model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected), r.RooFit.Components("topext"),\
    r.RooFit.LineColor(2), r.RooFit.LineStyle(r.kDashed), r.RooFit.FillColor(2))
model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected), r.RooFit.Components("wjetsext"),\
    r.RooFit.LineColor(r.kGreen-1), r.RooFit.LineStyle(r.kDashed))
model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected), r.RooFit.Components("wwext"),\
    r.RooFit.LineColor(r.kRed+1), r.RooFit.LineStyle(r.kDashed))
model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected), r.RooFit.Components("zzext"),\
    r.RooFit.LineColor(r.kYellow+1), r.RooFit.LineStyle(r.kDashed))
model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected), r.RooFit.Components("zext"),\
    r.RooFit.LineColor(r.kBlue+1), r.RooFit.LineStyle(r.kDashed))
f.Draw()
l1 = r.TLegend(0.1,0.7,0.4,0.9)
l1.AddEntry("data","Data", "lep")
e = l1.AddEntry("topext","t\\bar{t} + tW","l")
e.SetLineColor(2)
e.SetLineStyle(r.kDashed)
e.SetLineWidth(3)
e = l1.AddEntry("wjetsext","W+Jets","l")
e.SetLineColor(r.kGreen-1)
e.SetLineStyle(r.kDashed)
e.SetLineWidth(3)
e = l1.AddEntry("wwext","WW","l")
e.SetLineColor(r.kRed+1)
e.SetLineStyle(r.kDashed)
e.SetLineWidth(3)
e = l1.AddEntry("zzext","ZZ","l")
e.SetLineColor(r.kYellow+1)
e.SetLineStyle(r.kDashed)
e.SetLineWidth(3)
e=l1.AddEntry("zext","Z/#gamma^{*}","l")
e.SetLineColor(r.kBlue+1)
e.SetLineStyle(r.kDashed)
e.SetLineWidth(3)
l1.SetFillColor(r.kWhite)
l1.Draw()

# <markdowncell>

# Get the values out of the fit and compare to MC values

# <codecell>

ntop_low_mct = ntop.getVal()
ntop_low_mct_err = ntop.getError()
nwjets_low_mct = nwjets.getVal()
nwjets_low_mct_err = nwjets.getError()
nww_low_mct = nww.getVal()
nww_low_mct_err = nww.getError()
nzz_low_mct = nzz.getVal()
nzz_low_mct_err = nzz.getError()
nz_low_mct = nz.getVal()
nz_low_mct_err = nz.getError()
print ntop_low_mct_true, ntop_low_mct, ntop_low_mct_err
print nwjets_low_mct_true, nwjets_low_mct, nwjets_low_mct_err
print nww_low_mct_true, nww_low_mct, nww_low_mct_err
print nzz_low_mct_true, nzz_low_mct, nzz_low_mct_err
print nz_low_mct_true, nz_low_mct, nz_low_mct_err

# <markdowncell>

# Look at the results

# <codecell>

print "Actual Background: ",  sum(data[sel['sig_mct_high']].weight)
print "Predicted Total:", ntot_sig_pred.getVal(), "+-", ntot_sig_pred.getPropagatedError(results)
print "Actual Top: ",  sum(data[sel['sig_mct_high'] & ((data.mctype=="ttbar") | (data.mctype=='tW'))].weight)
print "Predicted Top:", ntop_sig_pred.getVal(), "+-", ntop_sig_pred.getPropagatedError(results)
print "Actual WW: ",  sum(data[sel['sig_mct_high'] & ((data.mctype=="WW") | (data.mctype=="WZ"))].weight)
print "Predicted WW:", nww_sig_pred.getVal(), "+-", nww_sig_pred.getPropagatedError(results)
print "Actual ZZ: ",  sum(data[sel['sig_mct_high'] & (data.mctype=="ZZ")].weight)
print "Predicted ZZ:", nzz_sig_pred.getVal(), "+-", nzz_sig_pred.getPropagatedError(results)
print "Actual Z: ",  sum(data[sel['sig_mct_high'] & (data.mctype=="DY")].weight)
print "Predicted Z:", nz_sig_pred.getVal(), "+-", nz_sig_pred.getPropagatedError(results)
print "Actual W: ",  sum(data[sel['sig_mct_high'] & (data.mctype=="wjets")].weight)
print "Predicted W:", nwjets_sig_pred.getVal(), "+-", nwjets_sig_pred.getPropagatedError(results)

# <codecell>


