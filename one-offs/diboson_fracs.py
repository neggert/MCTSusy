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

	total, bins  = np.histogram(vv.mctperp, weights=vv.weight, bins=19, range=(10,200))
	total[-1] += vv[vv.mctperp > 200].weight.sum()
	ww_hist, _ = np.histogram(ww.mctperp, weights=ww.weight, bins=19, range=(10,200))
	ww_hist[-1] += ww[ww.mctperp > 200].weight.sum()
	wz_hist, _ = np.histogram(wz.mctperp, weights=wz.weight, bins=19, range=(10,200))
	wz_hist[-1] += wz[wz.mctperp > 200].weight.sum()
	zz_hist, _ = np.histogram(zz.mctperp, weights=zz.weight, bins=19, range=(10,200))
	zz_hist[-1] += zz[zz.mctperp > 200].weight.sum()
	rare_hist, _ = np.histogram(rare.mctperp, weights=rare.weight, bins=19, range=(10,200))
	rare_hist[-1] += rare[rare.mctperp > 200].weight.sum()

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