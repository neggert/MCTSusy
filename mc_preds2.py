from config.data import *
from selection import *
from collections import defaultdict
import simplejson as json
import numpy as np

preds = defaultdict(dict)

chans = ['of', 'sf']

wz_fs = (mc.mc_cat=="WZ") & abs(mc.parentParentPdg1).isin([22,23]) & abs(mc.parentParentPdg1).isin([22,23])

mczz = mc[(mc.mc_cat=='ZZ')]
mcwz = mc[wz_fs]
mcvv = mczz.append(mcwz)
selvv = selection.get_samples( mcvv, 100.)

for ch in chans:
    preds[ch]['top'] = mc[smc['sig_mct_low_sf'] & (((mc.mc_cat=="top") & (mc.gen_neutrinos>=2) ) | (mc.mc_cat=="WW") | wz_fs | (mc.mc_cat=="HWW") | (mc.mc_cat=="VVV"))].weight.sum()

    preds[ch]['z'] = mc[smc['sig_mct_low_'+ch] & (mc['mc_cat']=='DY')].weight.sum()
    preds[ch]['vv'] = mc[smc['sig_mct_low_'+ch] & ((mc.mc_cat=="ZZ") | ((mc.mc_cat=="WZ") & ~(abs(mc.parentParentPdg1).isin([22,23]) & abs(mc.parentParentPdg1).isin([22,23]))))].weight.sum()
    preds[ch]['fake'] = mc[smc['sig_mct_low_'+ch] & ((mc['mc_cat']=='fake') | ((mc['mc_cat']=='top') & (mc['gen_neutrinos'] < 2)))].weight.sum()
    preds[ch]['sum'] = mc[smc['sig_mct_low_'+ch]].weight.sum()

print preds

with open("mc_preds2.json", 'w') as f:
    json.dump(preds, f)