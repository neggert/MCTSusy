#!/usr/bin/env python
"""
Combine hypothesis test results

Usage:
    find_limits.py [-p] <output_file> <input_files>... 

Options:
    -p          Make plots

"""

import sys
import re
import ROOT as R
import docopt
import scipy.interpolate
import scipy.optimize
import scipy.stats
import sklearn.gaussian_process as gp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

import warnings
warnings.simplefilter("ignore", DeprecationWarning)



def main(input_files, output_file, plot):

    tests = {}

    for f in input_files:
        m1, m2 = map(int, re.search("_(\d+)_(\d+)[_\.]", f).groups())
        rf = R.TFile(f)
        if (m1,m2) not in tests.keys():
            tests[(m1,m2)] = rf.Get("result_sig_strength")
        else:
            tests[(m1,m2)].Add(rf.Get("result_sig_strength"))

    with open(output_file, 'w') as f:
        for masses in tests.keys():
            do_plot = plot
            res = tests[masses]
            m1,m2 = masses

            pois = []
            obs_CLs_vs_poi = []
            exp_CLs_vs_poi = []
            expp1_CLs_vs_poi = []
            expm1_CLs_vs_poi = []
            for i in xrange(res.ArraySize()):
                poi = res.GetXValue(i)
                obs_cls = res.CLs(i)
                exp_cls_dist = res.GetExpectedPValueDist(i)
                size = exp_cls_dist.GetSize()
                # if size > 2000:
                #     print "{0} toys at point".format(size), poi, m1, m2
                exp_cls = exp_cls_dist.InverseCDF(0.5)
                expp1_cls = exp_cls_dist.InverseCDF(0.84)
                expm1_cls = exp_cls_dist.InverseCDF(0.16)
                if any(map(np.isnan, [obs_cls, exp_cls, expp1_cls, expm1_cls])):
                    raise ValueError("Found NaN CLs.")
                pois.append(poi)
                obs_CLs_vs_poi.append(obs_cls)
                exp_CLs_vs_poi.append(exp_cls)
                expp1_CLs_vs_poi.append(expp1_cls)
                expm1_CLs_vs_poi.append(expm1_cls)

            sort_order = np.argsort(pois)
            pois = np.asarray(pois)[:,np.newaxis][sort_order]
            obs_CLs_vs_poi = np.asarray(obs_CLs_vs_poi)[sort_order]
            exp_CLs_vs_poi = np.asarray(exp_CLs_vs_poi)[sort_order]
            expp1_CLs_vs_poi = np.asarray(expp1_CLs_vs_poi)[sort_order]
            expm1_CLs_vs_poi = np.asarray(expm1_CLs_vs_poi)[sort_order]
            target_cls = 0.05

            if max(pois) > 999 and any([min(obs_CLs_vs_poi)>target_cls, min(exp_CLs_vs_poi)>target_cls,
                                        min(expp1_CLs_vs_poi)>target_cls, min(expm1_CLs_vs_poi)>target_cls]):
                continue

            try: 
                obs_limit, obs, do_plot_a = get_limit(pois, obs_CLs_vs_poi, target_cls, m1, m2)
                exp_limit, exp, do_plot_b = get_limit(pois, exp_CLs_vs_poi, target_cls, m1, m2)
                expp1_limit, expp1, do_plot_c = get_limit(pois, expp1_CLs_vs_poi, target_cls, m1, m2)
                expm1_limit, expm1, do_plot_d = get_limit(pois, expm1_CLs_vs_poi, target_cls, m1, m2)
            except ValueError:
                print "Not enough non-zero points at ", m1, m2
                print zip(pois, obs_CLs_vs_poi)
                print zip(pois, exp_CLs_vs_poi)


                continue
            # if any([do_plot_a, do_plot_b, do_plot_c, do_plot_d]) or do_plot:

	    #     fig = plt.figure()
            #     plt.plot(pois, obs_CLs_vs_poi)
            #     plt.plot(pois, exp_CLs_vs_poi)
            #     plt.plot(pois, expp1_CLs_vs_poi)
            #     plt.plot(pois, expm1_CLs_vs_poi)
            #     plt.xlabel("Signal Strength")
            #     plt.ylabel("CLs")
            #     plt.yscale("log")
            #     plt.axhline(0.05)
            #     plt.legend(['Observed', 'Expected', "+1 sigma", "-1 sigma"])
            #     plt.savefig("figs/{0}_{1}_cls.png".format(m1, m2))
            #     fig = plt.figure()

            #     cls_all = np.linspace(min(pois), max(pois), 100)
            #     plt.plot(cls_all, map(obs_limit.interp.predict, cls_all))
            #     plt.plot(cls_all, map(exp_limit.interp.predict, cls_all))
            #     plt.plot(cls_all, map(expp1_limit.interp.predict, cls_all))
            #     plt.plot(cls_all, map(expm1_limit.interp.predict, cls_all))
            #     plt.xlabel("Signal Strength")
            #     plt.ylabel("CLs")
            #     plt.yscale("log")
            #     plt.axhline(0.05)
            #     plt.legend(['Observed', 'Expected', "+1 sigma", "-1 sigma"])
            #     plt.savefig("figs/{0}_{1}_cls_interpolated.png".format(m1, m2))

            f.write("{0}\t{1}\t{2:.3f}\t{3:.3f}\t{4:.3f}\t{5:.3f}\n".format(m1, m2, float(exp), float(expm1), float(expp1), float(obs)))

class LimitFinder(object):
    """docstring for LimitFinder"""
    def __init__(self, poi, cls, target_cls):
        # self.get_non_zero_cls(poi, cls)
        self.poi = poi
        self.cls = cls
        self.target_cls = target_cls
        self.interp = gp.GaussianProcess()
        self.interp.fit(self.poi, self.cls)

    def get_non_zero_cls(self, poi, cls):
        non_zero_indices = np.where(cls > 0.)
        # add the first 0 index
        try:
            non_zero_indices = np.append(non_zero_indices, np.min(np.where(cls == 0)))
        except ValueError:
            # probably weren't any 0 points
            pass
        if len(non_zero_indices) < 2:
            raise ValueError("Need at least two non-zero CLs values")
        self.poi = poi[non_zero_indices]
        self.cls = cls[non_zero_indices]

    def limit(self):
        return scipy.optimize.brentq(lambda x: self.interp.predict(x)-self.target_cls, min(self.poi), max(self.poi))

    def get_rel_uncertainty_at_limit(self):
        _, unc = self.interp.predict(self.limit(), eval_MSE=True)
        return unc/self.limit()
        
def get_limit(poi, cls, target_cls, m1, m2):
    do_plot = False
    limit_finder = LimitFinder(poi, cls, target_cls)
    try:
        limit = limit_finder.limit()
        if limit_finder.get_rel_uncertainty_at_limit() > 0.01:
            print "Large uncertainty at point", m1, m2
            print zip(poi, cls)
            do_plot = True
    except ValueError:
        print "Failed to interpolate limit for point", m1, m2
        print zip(poi, cls)
        do_plot=True
        limit = 0
        try:
            points = suggest_extra_points(poi, cls)
            for p in points:
                print get_file_name(m1, m2), p
        except:
            print "Couldn't suggest points"

    return limit_finder, limit, do_plot

def suggest_extra_points(poi, cls):
    if all(cls < 0.01):
        raise RuntimeWarning("all 0s")
    elif all(cls > 0.5):
        raise RuntimeWarning("Range is way off")
    elif min(cls) > 0.05:
        return extend_poi_higher(poi, cls)
    elif max(cls) < 0.05:
        return extend_poi_lower(poi, cls)
    elif not any((cls > 0.0) & (cls < 0.05)):
        return fill_poi(poi, cls)
    elif (np.diff(cls) > 0).any() :
        return rerun_bad_points(poi, cls)

def extend_poi_higher(poi, cls):
    if max(poi) == 1000.:
        raise RuntimeWarning("Out of Range")
    slope, intercept, r, p, err = scipy.stats.linregress(cls[-3:], poi[:,0][-3:])
    lower_end = max(poi)
    upper_end = (intercept-lower_end)*2+lower_end
    return np.linspace(lower_end, upper_end, 5)

def extend_poi_lower(poi, cls):
    slope, intercept, r, p, err = scipy.stats.linregress(cls[3:], poi[:,0][3:])
    upper_end = min(poi)
    lower_end = upper_end-(upper_end-intercept)*2
    return np.linspace(lower_end, upper_end, 5)

def fill_poi(poi, cls):
    sorted_ids = np.argsort(poi)
    s_poi = poi[sorted_ids]
    s_cls = cls[sorted_ids]
    indices = np.where((s_cls < 0.05) & (s_cls > 0.0))
    poi_bracket = s_poi[indices]
    poi_bracket = np.append(poi_bracket, s_poi[max(indices)+1])
    lower_end = min(poi_bracket)
    upper_end = max(poi_bracket)
    return np.linspace(lower_end, upper_end, 5)

def rerun_bad_points(poi, cls):
    ids = np.where(np.diff(cls) > 0)
    return poi[ids, 0].reshape(len(ids))

def get_file_name(m1, m2):
    return "../limits/sig_chi_{}_{}_combined_meas_model.root".format(m1,m2)

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    files = args['<input_files>']

    plot = bool(args['-p'])
    output_file = args['<output_file>']
    main(files, output_file, plot)

