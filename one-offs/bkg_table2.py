import json

from config.data import *

with open("fit_results2.json") as f:
	r = json.load(f)

dibosons = ['WW', 'WZ', 'ZZ', 'HWW','VVV']

vv_sf = mc[smc['sig_mct_low_sf'] & (mc.mc_cat.isin(dibosons))]

all_dibosons_sf = vv_sf.weight.sum()

WW_sf_frac = vv_sf[vv_sf.mc_cat=="WW"].weight.sum()/all_dibosons_sf
WZ_sf_frac = vv_sf[vv_sf.mc_cat=="WZ"].weight.sum()/all_dibosons_sf
ZZ_sf_frac = vv_sf[vv_sf.mc_cat=="ZZ"].weight.sum()/all_dibosons_sf
Rare_sf_frac = vv_sf[(vv_sf.mc_cat=="VVV") | (vv_sf.mc_cat=="HWW")].weight.sum()/all_dibosons_sf

print "sf"
print "WW", [x*WW_sf_frac for x in r['sf']['vv']]
print "WZ", [x*WZ_sf_frac for x in r['sf']['vv']]
print "ZZ", [x*ZZ_sf_frac for x in r['sf']['vv']]
print "Rare", [x*Rare_sf_frac for x in r['sf']['vv']]