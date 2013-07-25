import json

from config.data import *

with open("fit_results.json") as f:
	r = json.load(f)

dibosons = ['WW', 'WZ', 'ZZ', 'HWW','VVV']
vv_of = mc[smc['sig_mct_low_of'] & (mc.mc_cat.isin(dibosons))]

all_dibosons_of = vv_of.weight.sum()

WW_of_frac = vv_of[vv_of.mc_cat=="WW"].weight.sum()/all_dibosons_of
WZ_of_frac = vv_of[vv_of.mc_cat=="WZ"].weight.sum()/all_dibosons_of
ZZ_of_frac = vv_of[vv_of.mc_cat=="ZZ"].weight.sum()/all_dibosons_of
Rare_of_frac = vv_of[(vv_of.mc_cat=="VVV") | (vv_of.mc_cat=="HWW")].weight.sum()/all_dibosons_of

print "of"
print "WW", [x*WW_of_frac for x in r['of']['vv']]
print "WZ", [x*WZ_of_frac for x in r['of']['vv']]
print "ZZ", [x*ZZ_of_frac for x in r['of']['vv']]
print "Rare", [x*Rare_of_frac for x in r['of']['vv']]


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