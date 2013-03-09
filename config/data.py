from pandas import *
from config.parameters import *
import sys
import selection
reload(selection)
import numpy as np
import ROOT as R

# background MC
s = HDFStore("work/mc/mc_20130308.hdf5")
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

# pileup reweighting
mc_pu_hist, _ = np.histogram(mc.nTruePuVertices, bins=101, range=(0, 101))
mc_pu_hist = np.asarray(mc_pu_hist, dtype=np.float)

pu_file = R.TFile("config/TruePU.root")
pu_th1 = pu_file.Get("pileup")
data_pu_hist = np.zeros(mc_pu_hist.shape)

for i in xrange(len(data_pu_hist)):
	data_pu_hist[i] = pu_th1.GetBinContent(i)


# normalize the histograms
mc_pu_hist *= 1./np.sum(mc_pu_hist)
data_pu_hist *= 1./np.sum(data_pu_hist)

# calculate weights
pu_weights = data_pu_hist/mc_pu_hist

# apply the weights
mc.nTruePuVertices[mc.nTruePuVertices > 100] = 100
event_pu_weights = mc.nTruePuVertices.apply(lambda n: pu_weights.item(int(n)))
mc.weight *= event_pu_weights

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

sslep = HDFStore("work/sms/sms_slep.hdf5")
slep = sslep['data']
sel_slep = selection.get_samples(slep)
mumu_high_eta_slep = sel_slep['opposite_sign_mumu'] & (abs(slep.eta2) > 1.)
mumu_low_eta_slep = sel_slep['opposite_sign_mumu'] & (abs(slep.eta2) < 1.)

# slep.weight *= (sel_slep['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_slep.astype(float)*mumu_high_eta_trigger_eff
              # +mumu_low_eta_slep.astype(float)*mumu_low_eta_trigger_eff + sel_slep['opposite_sign_emu'].astype(float)*emu_trigger_eff)

stchiww = HDFStore("work/sms/sms_tchiww.hdf5")
tchiww = stchiww['data']
sel_tchiww = selection.get_samples(tchiww)
