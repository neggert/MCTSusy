from pandas import *
from config.parameters import *
import sys
import selection
reload(selection)
import numpy as np
import ROOT as R

# background MC
s = HDFStore("work/mc/mc_20131227.hdf5")
mc = s['data']
mc = mc[mc.mctype != "WGStarToLNu2E"]

weights = mc.x_eff*lumi
weights.name="weight"

mc = mc.join(weights)

# apply trigger efficiencies

smc = selection.get_samples(mc)

mumu_high_eta = smc['mumu'] & (abs(mc.eta2) > 1.)
mumu_low_eta = smc['mumu'] & (abs(mc.eta2) < 1.)


mc.weight *= (smc['ee'].astype(float)*ee_trigger_eff+mumu_high_eta.astype(float)*mumu_high_eta_trigger_eff
              +mumu_low_eta.astype(float)*mumu_low_eta_trigger_eff + smc['emu'].astype(float)*emu_trigger_eff)

# adjust some MC categories
# cat = mc.pop('mc_cat')
# cat[mc.mctype=="WWZNoGstar"] = "VVV"
# cat[mc.mctype=="WWW"] = 'VVV'
# mc = mc.join(cat)

# b-tag reweighting
jet0 = mc[['jet0pt', 'jet0_istagged']][~np.isnan(mc['jet0pt'])]
jet1 = mc[['jet1pt', 'jet1_istagged']][~np.isnan(mc['jet1pt'])]
jet2 = mc[['jet2pt', 'jet2_istagged']][~np.isnan(mc['jet2pt'])]
jet0.columns = ['jetpt', 'istagged']
jet1.columns = ['jetpt', 'istagged']
jet2.columns = ['jetpt', 'istagged']

jets = jet0.append(jet1, ignore_index=True).append(jet2, ignore_index=True)

def eff(df):
    return 1. * df[df.istagged].istagged.count() / df.istagged.count()

bins = np.arange(20, 810, 10)
effs = jets.groupby(np.digitize(jets.jetpt, bins)).agg(eff)['jetpt']

effs[effs > 0.99] = 0.99

def get_eff(pt):
    return effs[int(pt) / 10]

def sf(x):
    return (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x))

def get_weight(data):
    w = 1.
    for i in xrange(3):
        pt = data['jet'+str(i)+'pt']
        if not pt:
            continue
        if pt < 20 or pt > 800 or np.isnan(pt):
            continue
        if data['jet'+str(i)+'_istagged']:
            w *= sf(pt)
        else:
            w *= (1 - get_eff(pt)*sf(pt)) / (1 - get_eff(pt))

    if (np.isnan(w)):
        raise RuntimeError()
    return w

b_weights = mc.apply(get_weight, axis=1)

mc['weight'] *= b_weights

# pileup reweighting
mc_pu_hist = np.asarray([                         2.560E-06,
                         5.239E-06,
                         1.420E-05,
                         5.005E-05,
                         1.001E-04,
                         2.705E-04,
                         1.999E-03,
                         6.097E-03,
                         1.046E-02,
                         1.383E-02,
                         1.685E-02,
                         2.055E-02,
                         2.572E-02,
                         3.262E-02,
                         4.121E-02,
                         4.977E-02,
                         5.539E-02,
                         5.725E-02,
                         5.607E-02,
                         5.312E-02,
                         5.008E-02,
                         4.763E-02,
                         4.558E-02,
                         4.363E-02,
                         4.159E-02,
                         3.933E-02,
                         3.681E-02,
                         3.406E-02,
                         3.116E-02,
                         2.818E-02,
                         2.519E-02,
                         2.226E-02,
                         1.946E-02,
                         1.682E-02,
                         1.437E-02,
                         1.215E-02,
                         1.016E-02,
                         8.400E-03,
                         6.873E-03,
                         5.564E-03,
                         4.457E-03,
                         3.533E-03,
                         2.772E-03,
                         2.154E-03,
                         1.656E-03,
                         1.261E-03,
                         9.513E-04,
                         7.107E-04,
                         5.259E-04,
                         3.856E-04,
                         2.801E-04,
                         2.017E-04,
                         1.439E-04,
                         1.017E-04,
                         7.126E-05,
                         4.948E-05,
                         3.405E-05,
                         2.322E-05,
                         1.570E-05,
                         5.005E-06])

pu_file = R.TFile("config/TruePU.root")
pu_th1 = pu_file.Get("pileup")
data_pu_hist = np.zeros(mc_pu_hist.shape)

for i in xrange(len(data_pu_hist)):
	data_pu_hist[i] = pu_th1.GetBinContent(i+1)


# normalize the histograms
mc_pu_hist *= 1./np.sum(mc_pu_hist)
data_pu_hist *= 1./np.sum(data_pu_hist)

# calculate weights
pu_weights = data_pu_hist/mc_pu_hist

# apply the weights
mc.nTruePuVertices[mc.nTruePuVertices > 60] = 60
event_pu_weights = mc.nTruePuVertices.apply(lambda n: pu_weights.item(int(n)))
mc.weight *= event_pu_weights

# Z MC MET re-weighting
import json
with open("z_weights.json") as f:
	reweighting = json.load(f)

z_weights = np.asarray(reweighting['scale_factors'])
z_bins = np.asarray(reweighting['bins'])

def reweight(val, weight_factors, bins):
	i = np.argmax(np.where(bins<val, bins, 0))
	if i > len(weight_factors)-1:
		return 1.
	return weight_factors[i]

mc.weight[mc.mc_cat=="DY"] *= mc[mc.mc_cat=="DY"].metPt.apply(reweight, args=(z_weights, z_bins))

# data
t = HDFStore("work/data/data_20130304.hdf5")
data = t['data']
sd = selection.get_samples(data, 100., True)
# add a weight column
data = data.join(Series(np.ones(data.mctperp.count()), name="weight", index=data.index))

# signal MC
schi = HDFStore("work/sms/sms_chi_25GeV.hdf5")
chi = schi['data']
sel_chi = selection.get_samples(chi)
mumu_high_eta_chi = sel_chi['opposite_sign_mumu'] & (abs(chi.eta2) > 1.)
mumu_low_eta_chi = sel_chi['opposite_sign_mumu'] & (abs(chi.eta2) < 1.)

weight = (sel_chi['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_chi.astype(float)*mumu_high_eta_trigger_eff
              +mumu_low_eta_chi.astype(float)*mumu_low_eta_trigger_eff + sel_chi['opposite_sign_emu'].astype(float)*emu_trigger_eff)
weight.name="weight"
chi = chi.join(weight)

# FastSim -> FullSim reweighting
mu_weight_file = R.TFile("config/muon_FastSim_EWKino.root")
mu_fastsim_weights = mu_weight_file.Get("SF")
ele_weight_file = R.TFile("config/electron_FastSim_EWKino.root")
ele_fastsim_weights = mu_weight_file.Get("SF")

def fastsim_weight(row):
    if abs(row['pdg1']) == 13:
        s1 = mu_fastsim_weights.GetBinContent(mu_fastsim_weights.GetXaxis().FindBin(row['pt1']), mu_fastsim_weights.GetYaxis().FindBin(abs(row['eta1'])))
    else:
        s1 = ele_fastsim_weights.GetBinContent(ele_fastsim_weights.GetXaxis().FindBin(row['pt1']), ele_fastsim_weights.GetYaxis().FindBin(abs(row['eta1'])))
    if abs(row['pdg2']) == 13:
        s2 = mu_fastsim_weights.GetBinContent(mu_fastsim_weights.GetXaxis().FindBin(row['pt2']), mu_fastsim_weights.GetYaxis().FindBin(abs(row['eta2'])))
    else:
        s2 = ele_fastsim_weights.GetBinContent(ele_fastsim_weights.GetXaxis().FindBin(row['pt2']), ele_fastsim_weights.GetYaxis().FindBin(abs(row['eta2'])))

    return s1 * s2

chi.weight *= chi.apply(fastsim_weight, axis=1)

# add PU reweighting

sslep = HDFStore("work/sms/sms_slep.hdf5")
slep = sslep['data']
sel_slep = selection.get_samples(slep)
mumu_high_eta_slep = sel_slep['opposite_sign_mumu'] & (abs(slep.eta2) > 1.)
mumu_low_eta_slep = sel_slep['opposite_sign_mumu'] & (abs(slep.eta2) < 1.)

weight = (sel_slep['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_slep.astype(float)*mumu_high_eta_trigger_eff
              +mumu_low_eta_slep.astype(float)*mumu_low_eta_trigger_eff + sel_slep['opposite_sign_emu'].astype(float)*emu_trigger_eff)
weight.name="weight"

slep = slep.join(weight)

slep.weight *= slep.apply(fastsim_weight, axis=1)


stchiww = HDFStore("work/sms/sms_tchiww.hdf5")
tchiww = stchiww['data']
sel_tchiww = selection.get_samples(tchiww)
