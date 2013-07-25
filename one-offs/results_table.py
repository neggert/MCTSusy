from prettytable import PrettyTable

from config.data import *

import json

with open(sys.argv[1]) as f:
	results = json.load(f)

t = PrettyTable()
t.vertical_char = "&"

vv_cats = ['WW', 'WZ', 'ZZ', 'HWW', 'VVV']
categories = ['top', 'WW', 'WZ', 'ZZ', 'HWW', 'VVV', 'DY', 'fake']
t.add_column("", categories)

for channel in ['of', 'sf']:
	sig_mc = mc[smc['sig_mct_low_'+channel]]
	vv_mc = sig_mc[sig_mc.mc_cat.isin(vv_cats)]
	cola = []
	colb = []
	for bkg in categories:
		if bkg in vv_cats:
			frac = sig_mc[sig_mc.mc_cat==bkg].weight.sum()/vv_mc.weight.sum()
			count, err = results[channel]['vv']
			count *= frac
			err *= frac
		else:
			count, err = results[channel][bkg]

		mc_count = float(sig_mc[sig_mc.mc_cat==bkg].weight.sum())
		if count > 100 or mc_count > 100 or err > 100:
			cola.append("${0:.0f}\pm{1:.0f}$".format(count, err))
			colb.append("{0:.0f}".format(mc_count))
		else:
			cola.append("${0:.2g}\pm{1:.2g}$".format(count, err))
			colb.append("{0:.2g}".format(mc_count))

	t.add_column(channel+"_fit", cola)
	t.add_column(channel+"_mc", colb)

print t



