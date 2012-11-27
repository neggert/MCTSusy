from pandas import *
from parameters import *
import sys
sys.path.append("../")
import selection
reload(selection)

s = HDFStore("Data/mc2012_201211108.hdf5")
mc = s['all_cat']
t = HDFStore("Data/good_data_triggers.hdf5")
data = t['data']

# apply trigger efficiencies

smc = selection.get_samples(mc)

mumu_high_eta = smc['opposite_sign_mumu'] & (abs(mc.eta2) > 1.)
mumu_low_eta = smc['opposite_sign_mumu'] & (abs(mc.eta2) < 1.)


mc.weight *= (smc['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta.astype(float)*mumu_high_eta_trigger_eff
              +mumu_low_eta.astype(float)*mumu_low_eta_trigger_eff + smc['opposite_sign_emu'].astype(float)*emu_trigger_eff)

sd = selection.get_samples(data, 100., False)

schi = HDFStore("Data/sms_chi.hdf5")
chi = schi['data']

sslep = HDFStore("Data/sms_slep.hdf5")
slep = sslep['data']
