from config.data import *

import numpy as np
import scipy.optimize as sopt

# z_range = (10, 60)
# z_bins = 5

# onz_data = data[sd['z_ctrl_sf']]
# onz_mc = mc[smc['z_ctrl_sf']]
# z_mc = onz_mc[onz_mc.mc_cat=="DY"]
# other_mc = onz_mc[onz_mc.mc_cat!="DY"]

# data_hist, bins = np.histogram(onz_data.mctperp, bins=z_bins, range=z_range)
# z_hist, bins = np.histogram(z_mc.mctperp, weights=z_mc.weight, bins=z_bins, range=z_range)
# other_hist, bins = np.histogram(other_mc.mctperp, weights=other_mc.weight, bins=z_bins, range=z_range)

# def chi2(fitted, data):
# 	return np.sum((fitted - data)**2/np.sqrt(data))

# def fit_func(alpha, z, other):
# 	return alpha*z+other

# def min_func(alpha, fit, other, data):
# 	return chi2(fit_func(alpha, fit, other), data)

# res = sopt.minimize_scalar(min_func, args=(z_hist, other_hist, data_hist))

# print res

onz_data = data[sd['z_ctrl_0met']]
onz_mc = mc[smc['z_ctrl_0met']]
z_mc = onz_mc[onz_mc.mc_cat=="DY"]
other_mc = onz_mc[onz_mc.mc_cat!="DY"]


met_bins = np.arange(0, z_mc.metPt.max(), 10)
met_bins = np.append(met_bins, z_mc.metPt.max())


finished = False

while not finished:
	data_hist, bins = np.histogram(onz_data.metPt, bins=met_bins)
	z_hist, bins = np.histogram(z_mc.metPt, weights=z_mc.weight, bins=met_bins)
	other_hist, bins = np.histogram(other_mc.metPt, weights=other_mc.weight, bins=met_bins)

	subbed_data_hist = data_hist-other_hist
	# subbed_data_hist[subbed_data_hist < 0 ] = 0.

	print subbed_data_hist

	# remove bin edges where we get < 5 counts in data
	if np.any(subbed_data_hist < 5):
		i = np.argmax(np.where(subbed_data_hist < 5, bins[:-1], 0))
		met_bins = np.delete(met_bins, i)
	elif np.any(z_hist < 5):
		i = np.argmax(np.where(z_hist < 5, bins[:-1], 0))
		met_bins = np.delete(met_bins, i)
	else:
		finished = True

print bins
# z_hist[z_hist<0.1] = 0

z_scales = subbed_data_hist/z_hist

import json

output = {"scale_factors":z_scales.tolist(),
		  "bins": bins.tolist()
		 }

with open("z_weights.json", 'w') as f:
	json.dump(output, f)

def reweight(val, weight_factors, bins):
	i = np.argmax(np.where(bins<val, bins, 0))
	if i > len(weight_factors)-1:
		return 1.
	return weight_factors[i]

# mc.weight[mc.mc_cat=="DY"] *= mc[mc.mc_cat=="DY"].metPt.apply(reweight, args=(z_scales, met_bins))