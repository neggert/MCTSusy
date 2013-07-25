from config.data import *
import numpy as np
import simplejson

chs = ['of', 'sf']

def get_fracs_ch( ch ):
	vv = mc[smc['sig_mct_low_'+ch] & mc.mc_cat.isin(['WW','WZ', 'ZZ', 'VVV', 'HWW'])]
	ww = vv[vv.mc_cat=='WW']
	wz = vv[vv.mc_cat=='WZ']
	zz = vv[vv.mc_cat=='ZZ']
	rare = vv[vv.mc_cat.isin(['VVV', 'HWW'])]

	total, bins  = np.histogram(vv.mctperp, weights=vv.weight, bins=29, range=(10,300))
	ww_hist, _ = np.histogram(ww.mctperp, weights=ww.weight, bins=29, range=(10,300))
	wz_hist, _ = np.histogram(wz.mctperp, weights=wz.weight, bins=29, range=(10,300))
	zz_hist, _ = np.histogram(zz.mctperp, weights=zz.weight, bins=29, range=(10,300))
	rare_hist, _ = np.histogram(rare.mctperp, weights=rare.weight, bins=29, range=(10,300))

	out = {'WW': (1.*ww_hist/total).tolist(),
		   'WZ': (1.*wz_hist/total).tolist(),
		   'ZZ': (1.*zz_hist/total).tolist(),
		   'Rare': (1.*rare_hist/total).tolist()
		  }
	return out

if __name__ == '__main__':
	out = {}
	for ch in chs:
		out[ch] = get_fracs_ch(ch)

	with open("diboson_fracs.json", 'w') as f:
		json.dump(out, f)