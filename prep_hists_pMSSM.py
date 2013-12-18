#!/usr/bin/env python

"""Prepare histograms

Usage:
    prep_hists.py <signal_file> <signal_output> [-b bins] [-l LOW] [-u high] [-x xsec]

Options:
    -h --help               Show this help information
    -b --bins=<b>              Number of bins [default: 19]
    -l --low=<l>               Lower end of histograms [default: 10]
    -u --high=<h>              Upper end of histograms [default: 200]
    -x --xsec_multiplier=<x>   Cross-section multiplier [default: 1]
"""

from background_prediction import rootutils
reload(rootutils)
import selection
reload(selection)

from config.parameters import lumi
from config.data import *

import sys
sys.path.append("import")
import CMSPyLibs.general_calc as general_calc

import ROOT as R
import os
import numpy as np

channels = ['of', 'sf']
backgrounds = ['of', 'vv', 'wjets', 'z']

import json

def load_xsec(filename):
    """
    load cross-sections for a file. Returns a dictionary describing cross-sections on their uncertainties.
    The dictionary is indexed by the particle mass, and each element is a tuple containing (xsec, uncertainty)
    """
    rf = R.TFile("limits/pMSSMtree_lhd.root")
    t = rf.Get("pMSSM")
    xsec_dict = {}
    for row in t:
        xsec_dict[row.ID] = row.xsect_pb
    rf.Close()
    return xsec_dict


def create_signal_file(input_file, out_filename, hist_filename, xsec_filename, xsec_multiplier=1., bins=19, histrange=(10,200) ):

    sms_file = HDFStore(input_file)
    sms = sms_file['data']
    sel_sms = selection.get_samples(sms)
    mumu_high_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) > 1.)
    mumu_low_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) < 1.)

    weight = (sel_sms['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_sms.astype(float)*mumu_high_eta_trigger_eff
                  +mumu_low_eta_sms.astype(float)*mumu_low_eta_trigger_eff + sel_sms['opposite_sign_emu'].astype(float)*emu_trigger_eff)
    weight.name="weight"
    sms = sms.join(weight)

    # pu reweighting
    # sms_pu_hist, _ = np.histogram(sms.nTruePuVertices, bins=101, range=(0, 101))
    # sms_pu_hist = np.asarray(sms_pu_hist, dtype=np.float)

    sms_pu_hist = np.asarray([2.344E-05,
                                  2.344E-05,
                                  2.344E-05,
                                  2.344E-05,
                                  4.687E-04,
                                  4.687E-04,
                                  7.032E-04,
                                  9.414E-04,
                                  1.234E-03,
                                  1.603E-03,
                                  2.464E-03,
                                  3.250E-03,
                                  5.021E-03,
                                  6.644E-03,
                                  8.502E-03,
                                  1.121E-02,
                                  1.518E-02,
                                  2.033E-02,
                                  2.608E-02,
                                  3.171E-02,
                                  3.667E-02,
                                  4.060E-02,
                                  4.338E-02,
                                  4.520E-02,
                                  4.641E-02,
                                  4.735E-02,
                                  4.816E-02,
                                  4.881E-02,
                                  4.917E-02,
                                  4.909E-02,
                                  4.842E-02,
                                  4.707E-02,
                                  4.501E-02,
                                  4.228E-02,
                                  3.896E-02,
                                  3.521E-02,
                                  3.118E-02,
                                  2.702E-02,
                                  2.287E-02,
                                  1.885E-02,
                                  1.508E-02,
                                  1.166E-02,
                                  8.673E-03,
                                  6.190E-03,
                                  4.222E-03,
                                  2.746E-03,
                                  1.698E-03,
                                  9.971E-04,
                                  5.549E-04,
                                  2.924E-04,
                                  1.457E-04,
                                  6.864E-05,
                                  3.054E-05,
                                  1.282E-05,
                                  5.081E-06,
                                  1.898E-06,
                                  6.688E-07,
                                  2.221E-07,
                                  6.947E-08,
                                  2.047E-08])

    pu_file = R.TFile("config/TruePU.root")
    pu_th1 = pu_file.Get("pileup")
    data_pu_hist = np.zeros(sms_pu_hist.shape)

    for i in xrange(len(data_pu_hist)):
        data_pu_hist[i] = pu_th1.GetBinContent(i+1)

    # normalize the histograms
    sms_pu_hist *= 1./np.sum(sms_pu_hist)
    data_pu_hist *= 1./np.sum(data_pu_hist)

    # calculate weights
    pu_weights = data_pu_hist/sms_pu_hist

    # apply the weights
    sms.nTruePuVertices[sms.nTruePuVertices > 60] = 60
    event_pu_weights = sms.nTruePuVertices.apply(lambda n: pu_weights.item(int(n)))
    sms.weight *= event_pu_weights

    # FastSim -> FullSim reweighting
    mu_weight_file = R.TFile("config/muon_FastSim_EWKino.root")
    mu_fastsim_weights = mu_weight_file.Get("SF")
    ele_weight_file = R.TFile("config/electron_FastSim_EWKino.root")
    ele_fastsim_weights = mu_weight_file.Get("SF")

    def fastsim_weight(row):
        if abs(row['pdg1']) == 13:
            s1 = mu_fastsim_weights.GetBinContent(mu_fastsim_weights.GetXaxis().FindBin(row['pt1']), mu_fastsim_weights.GetYaxis().FindBin(abs(row['eta1'])))
        else:
            s1 = ele_fastsim_weights.GetBinContent(ele_fastsim_weights.GetXaxis().FindBin(row['pt1']), ele_fastsim_weights.GetYaxis().FindBin(abs(row['eta1'])))
        if abs(row['pdg2']) == 13:
            s2 = mu_fastsim_weights.GetBinContent(mu_fastsim_weights.GetXaxis().FindBin(row['pt2']), mu_fastsim_weights.GetYaxis().FindBin(abs(row['eta2'])))
        else:
            s2 = ele_fastsim_weights.GetBinContent(ele_fastsim_weights.GetXaxis().FindBin(row['pt2']), ele_fastsim_weights.GetYaxis().FindBin(abs(row['eta2'])))

        return s1 * s2

    sms.weight *= sms.apply(fastsim_weight, axis=1)

    with open("limits/nevents_pMSSM.json") as f:
        nevents_dict = json.load(f)

    xsec_dict = load_xsec(xsec_filename)

    rfile = R.TFile(out_filename, "RECREATE")

    templates = {}
    jes_up_templates = {}
    jes_down_templates = {}

    masses = []

    # create a different template for each mass point
    # templates are scaled to number of exected events given the reference cross-section
    groups = sms.groupby(['mass1', 'mass2'])
    for name, data in groups:
        m1, m2 = name
        masses.append((m1, m2))
        sel = selection.get_samples(data)

        try:
            events_per_point = nevents_dict[str(m1)]
        except KeyError:
            print "Could not find events per point", m1
            continue
        try:
            xsec= xsec_dict[m1]
        except KeyError:
            print "Could not find xsection", m1

        xsec *= xsec_multiplier

        for ch in channels:
            mass_point = data[sel['sig_'+ch]]
            mass_point_jes_up = data[sel['sig_scaleup_'+ch]]
            mass_point_jes_down = data[sel['sig_scaledown_'+ch]]

            if mass_point[(mass_point.mctperp > histrange[0])].mctperp.count() == 0: continue

            h = rootutils.create_TH1(mass_point.mctperp, mass_point.weight,
                                                                         "sms_template_{}_{}_{}".format(ch, int(m1), int(m2)),
                                                                          bins, histrange, False)
            hup = rootutils.create_TH1(mass_point_jes_up.mctperp_up, mass_point_jes_up.weight,
                                                                         "sms_template_jes_up_{}_{}_{}".format(ch, int(m1), int(m2)),
                                                                          bins, histrange, False)
            hdown = rootutils.create_TH1(mass_point_jes_down.mctperp_down, mass_point_jes_down.weight,
                                                                         "sms_template_jes_down_{}_{}_{}".format(ch, int(m1), int(m2)),
                                                                          bins, histrange, False)
            h.Scale(xsec*lumi/events_per_point)
            hup.Scale(xsec*lumi/events_per_point)
            hdown.Scale(xsec*lumi/events_per_point)
            tname = 'sms_{}_{}_{}'.format(ch, m1, m2)
            templates[tname] = h
            jes_up_templates[tname] = hup
            jes_down_templates[tname] = hdown

    with open("pMSSM_points.json", 'w') as f:
        json.dump(masses, f)



    for k in templates.keys():
        templates[k].Write()
        jes_down_templates[k].Write()
        jes_up_templates[k].Write()

    rfile.Close()




if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)

    print(args)

    bins = int(args['--bins'])
    histrange = (float(args['--low']), float(args['--high']))

    create_signal_file(args['<signal_file>'], args['<signal_output>'], bins, histrange)

