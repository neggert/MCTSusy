#!/usr/bin/env python
"""
Combine hypothesis test results

Usage:
    combine_hypotest.py <output_file> <input_files>...

"""

import sys
import re
import ROOT as R
import docopt
import scipy.interpolate
import scipy.optimize
import sklearn.gaussian_process as gp
import numpy as np
import matplotlib.pylab as plt
plt.switch_backend("agg")

args = docopt.docopt(__doc__)

files = args['<input_files>']

tests = {}

for f in files:
    m1, m2 = map(int, re.search("_(\d+)_(\d+)_", f).groups())
    rf = R.TFile(f)
    if (m1,m2) not in tests.keys():
        tests[(m1,m2)] = rf.Get("result_sig_strength")
    else:
        tests[(m1,m2)].Add(rf.Get("result_sig_strength"))

with open(args['<output_file>'], 'w') as f:
    for masses in tests.keys():
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

        pois = np.asarray(pois)[:,np.newaxis]
        obs_CLs_vs_poi = np.asarray(obs_CLs_vs_poi)
        exp_CLs_vs_poi = np.asarray(exp_CLs_vs_poi)
        expp1_CLs_vs_poi = np.asarray(expp1_CLs_vs_poi)
        expm1_CLs_vs_poi = np.asarray(expm1_CLs_vs_poi)
        target_cls = 0.05

        try:
            interp_obs = gp.GaussianProcess()
            interp_obs.fit(pois, obs_CLs_vs_poi)
            obs = scipy.optimize.brentq(lambda x: interp_obs.predict(x)-target_cls, min(pois), max(pois))

            interp_exp = gp.GaussianProcess()
            interp_exp.fit(pois, exp_CLs_vs_poi)
            exp = scipy.optimize.brentq(lambda x: interp_exp.predict(x)-target_cls, min(pois), max(pois))

            interp_expp1 = gp.GaussianProcess()
            interp_expp1.fit(pois, expp1_CLs_vs_poi)
            expp1 = scipy.optimize.brentq(lambda x: interp_expp1.predict(x)-target_cls, min(pois), max(pois))

            interp_expm1 = gp.GaussianProcess()
            interp_expm1.fit(pois, expm1_CLs_vs_poi)
            expm1 = scipy.optimize.brentq(lambda x: interp_expm1.predict(x)-target_cls, min(pois), max(pois))
        except ValueError:
            continue


        # try:
        #     interp_obs = scipy.interpolate.UnivariateSpline(pois, obs_CLs_vs_poi)
        #     obs = scipy.optimize.brentq(lambda x: interp_obs(x)-target_cls, min(pois), max(pois))
        #     interp_exp = scipy.interpolate.UnivariateSpline(pois, exp_CLs_vs_poi)
        #     # import pdb; pdb.set_trace()
        #     exp = scipy.optimize.brentq(lambda x: interp_exp(x)-target_cls, min(pois), max(pois))
        #     interp_expp1 = scipy.interpolate.UnivariateSpline(pois, expp1_CLs_vs_poi)
        #     expp1 = scipy.optimize.brentq(lambda x: interp_expp1(x)-target_cls, min(pois), max(pois))
        #     interp_expm1 = scipy.interpolate.UnivariateSpline(pois, expm1_CLs_vs_poi)
        #     expm1 = scipy.optimize.brentq(lambda x: interp_expm1(x)-target_cls, min(pois), max(pois))
        # except ValueError:
        #     continue

        if any(map(np.isnan, [obs, exp, expp1, expm1])):
            fig = plt.figure()
            plt.plot(pois, obs_CLs_vs_poi)
            plt.plot(pois, exp_CLs_vs_poi)
            plt.plot(pois, expp1_CLs_vs_poi)
            plt.plot(pois, expm1_CLs_vs_poi)
            plt.xlabel("Signal Strength")
            plt.ylabel("CLs")
            plt.yscale("log")
            plt.axhline(0.05)
            plt.legend(['Observed', 'Expected', "+1 sigma", "-1 sigma"])
            plt.savefig("figs/{0}_{1}_cls.png".format(m1, m2))
            fig = plt.figure()

            cls_all = np.linspace(min(exp_CLs_vs_poi), max(exp_CLs_vs_poi), 100)
            plt.plot(map(interp_obs, cls_all), cls_all)
            plt.plot(map(interp_exp, cls_all), cls_all)
            plt.plot(map(interp_expp1, cls_all), cls_all)
            plt.plot(map(interp_expm1, cls_all), cls_all)
            plt.xlabel("Signal Strength")
            plt.ylabel("CLs")
            # plt.yscale("log")
            plt.axhline(0.05)
            plt.legend(['Observed', 'Expected', "+1 sigma", "-1 sigma"])
            plt.savefig("figs/{0}_{1}_cls_interpolated.png".format(m1, m2))


        # res.SetInterpolationOption(R.RooStats.HypoTestInverterResult.kSpline)
        # exp = res.GetExpectedUpperLimit()
        # exp_down = res.GetExpectedUpperLimit(-1)
        # exp_up = res.GetExpectedUpperLimit(+1)
        # obs = res.UpperLimit()


        f.write("{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4:.2f}\t{5:.2f}\n".format(m1, m2, float(exp), float(expm1), float(expp1), float(obs)))



