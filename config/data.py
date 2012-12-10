from pandas import *
from config.parameters import *
import sys
import selection
reload(selection)

s = HDFStore("work/mc/mc2012_201211108.hdf5")
mc = s['all_cat']
t = HDFStore("work/data/good_data_triggers.hdf5")
data = t['data']

# apply trigger efficiencies

smc = selection.get_samples(mc)

mumu_high_eta = smc['opposite_sign_mumu'] & (abs(mc.eta2) > 1.)
mumu_low_eta = smc['opposite_sign_mumu'] & (abs(mc.eta2) < 1.)


mc.weight *= (smc['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta.astype(float)*mumu_high_eta_trigger_eff
              +mumu_low_eta.astype(float)*mumu_low_eta_trigger_eff + smc['opposite_sign_emu'].astype(float)*emu_trigger_eff)

sd = selection.get_samples(data, 100., False)

schi = HDFStore("work/sms/sms_chi.hdf5")
chi = schi['data']
sel_chi = selection.get_samples(chi)
mumu_high_eta_chi = sel_chi['opposite_sign_mumu'] & (abs(chi.eta2) > 1.)
mumu_low_eta_chi = sel_chi['opposite_sign_mumu'] & (abs(chi.eta2) < 1.)

chi.weight *= (sel_chi['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_chi.astype(float)*mumu_high_eta_trigger_eff
              +mumu_low_eta_chi.astype(float)*mumu_low_eta_trigger_eff + sel_chi['opposite_sign_emu'].astype(float)*emu_trigger_eff)

sslep = HDFStore("work/sms/sms_slep.hdf5")
slep = sslep['data']
sel_slep = selection.get_samples(slep)
mumu_high_eta_slep = sel_slep['opposite_sign_mumu'] & (abs(slep.eta2) > 1.)
mumu_low_eta_slep = sel_slep['opposite_sign_mumu'] & (abs(slep.eta2) < 1.)

slep.weight *= (sel_slep['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_slep.astype(float)*mumu_high_eta_trigger_eff
              +mumu_low_eta_slep.astype(float)*mumu_low_eta_trigger_eff + sel_slep['opposite_sign_emu'].astype(float)*emu_trigger_eff)
