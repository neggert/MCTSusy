from config.data import *
from collections import defaultdict

import json

outdict = defaultdict(dict)

mc_cats = ['WW', 'WZ', 'ZZ', 'VVV', 'HWW']

for ch in ['of', 'sf']:
	vv = mc[smc['sig_mct_low_'+ch]&mc.mc_cat.isin(mc_cats)]
	for c in mc_cats:
		outdict[ch][c] = vv[vv.mc_cat==c].weight.sum()/vv.weight.sum()

with open("vv_fracs.json", 'w') as f:
	json.dump(outdict, f, indent=3)
	