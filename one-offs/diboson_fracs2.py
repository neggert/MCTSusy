from config.data import *
import numpy as np
import simplejson

chs = ['sf',]

def get_fracs_ch( ch ):
	vv = mc[((mc.mc_cat=="WZ") & abs(mc.parentParentPdg1).isin([22,23]) & abs(mc.parentParentPdg1).isin([22,23])) | (mc.mc_cat=="ZZ")]
	wz = vv[vv.mc_cat=='WZ']
	zz = vv[vv.mc_cat=='ZZ']

	total, bins  = np.histogram(vv.mctperp, weights=vv.weight, bins=29, range=(10,300))
	wz_hist, _ = np.histogram(wz.mctperp, weights=wz.weight, bins=29, range=(10,300))
	zz_hist, _ = np.histogram(zz.mctperp, weights=zz.weight, bins=29, range=(10,300))


	out = {
		   'WZ': (1.*wz_hist/total).tolist(),
		   'ZZ': (1.*zz_hist/total).tolist()
		  }
	return out

if __name__ == '__main__':
	out = {}
	for ch in chs:
		out[ch] = get_fracs_ch(ch)

	with open("diboson_fracs2.json", 'w') as f:
		json.dump(out, f)