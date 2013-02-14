from pandas import *
from config.parameters import *
import sys
import selection
reload(selection)

# background MC
s = HDFStore("work/mc/mc2012_moriond.hdf5")
mc = s['data']

weights = mc.x_eff*lumi
weights.name="weight"

mc = mc.join(weights)

# apply trigger efficiencies

smc = selection.get_samples(mc)

mumu_high_eta = smc['opposite_sign_mumu'] & (abs(mc.eta2) > 1.)
mumu_low_eta = smc['opposite_sign_mumu'] & (abs(mc.eta2) < 1.)


mc.weight *= (smc['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta.astype(float)*mumu_high_eta_trigger_eff
              +mumu_low_eta.astype(float)*mumu_low_eta_trigger_eff + smc['opposite_sign_emu'].astype(float)*emu_trigger_eff)

# adjust some MC categories
cat = mc.pop('mc_cat')
cat[mc.mctype=="WWZNoGstar"] = "VVV"
cat[mc.mctype=="WWW"] = 'VVV'
mc = mc.join(cat)


# data
t = HDFStore("work/data/data_final2012json.hdf5")
data = t['data']
sd = selection.get_samples(data, 100., False)

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
