from selection import *
from pandas import *
from config.data import *

def run_rare(filename):
	s = HDFStore(filename)
	d = s['data']
	sel = get_samples(d)

	mumu_high_eta = sel['opposite_sign_mumu'] & (abs(d.eta2) > 1.)
	mumu_low_eta = sel['opposite_sign_mumu'] & (abs(d.eta2) < 1.)

	weights = d.x_eff*lumi

	weights *= (sel['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta.astype(float)*mumu_high_eta_trigger_eff
	              +mumu_low_eta.astype(float)*mumu_low_eta_trigger_eff + sel['opposite_sign_emu'].astype(float)*emu_trigger_eff)

	weights.name="weight"
	d = d.join(weights)

	selected = d[sel['sig_mct_low']]

	return selected

