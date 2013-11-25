from config.data import *
from selection import *
from collections import defaultdict
import simplejson as json

preds = defaultdict(dict)

chans = ['of', 'sf']

for ch in chans:
    preds[ch]['top'] = mc[smc['sig_mct_low_'+ch] & (mc['mc_cat']=='top') & (mc['gen_neutrinos'] >= 2)].weight.sum()
    preds[ch]['z'] = mc[smc['sig_mct_low_'+ch] & (mc['mc_cat']=='DY')].weight.sum()
    preds[ch]['vv'] = mc[smc['sig_mct_low_'+ch] & (mc['mc_cat'].isin(['WW','ZZ','WZ','VVV','HWW']))].weight.sum()
    preds[ch]['fake'] = mc[smc['sig_mct_low_'+ch] & ((mc['mc_cat']=='fake') | ((mc['mc_cat']=='top') & (mc['gen_neutrinos'] < 2)))].weight.sum()
    preds[ch]['sum'] = mc[smc['sig_mct_low_'+ch]].weight.sum()

print preds

with open("mc_preds.json", 'w') as f:
    json.dump(preds, f)