#!/usr/bin/env python

"""Prepare histograms

Usage:
    prep_hists.py <signal_file> <signal_output> <xsec_file> <nevents_file> [-b bins] [-l LOW] [-u high] [-x xsec]

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


import ROOT as R
import os
import numpy as np

channels = ['of', 'sf']
backgrounds = ['of', 'vv', 'wjets', 'z']

def load_xsec(filename):
    """
    load cross-sections for a file. Returns a dictionary describing cross-sections on their uncertainties.
    The dictionary is indexed by the particle mass, and each element is a tuple containing (xsec, uncertainty)
    """
    f = open(filename)
    xsec_dict = {}
    for line in f:
        mass, xsec, err = line.split()
        xsec_dict[float(mass)] = (float(xsec), float(err))
    f.close()
    return xsec_dict

def create_template_file(filename="templates.root", bins=19, histrange=(10, 200)):
    """
    Create a ROOT file containing all of the background templates
    """

    mcvv = mc[(mc.mc_cat=='WW') | (mc.mc_cat=='ZZ') | (mc.mc_cat=='WZ')]
    mcz = mc[mc.mc_cat=='DY']
    selvv = selection.get_samples( mcvv, 100.)
    selz = selection.get_samples( mcz, 100.)

    rfile = R.TFile(filename, "RECREATE")
    rfile.cd()

    templates = {}

    constraints = {}

    for ch in channels:
        top = data[sd['top_ctrl_'+ch]]
        templates['top_'+ch] = rootutils.create_TH1(top.mctperp, top.weight, "top_template_"+ch, bins, histrange, True)

        wjets = data[sd['wjets_ctrl_'+ch]]
        templates['wjets_'+ch] = rootutils.create_TH1(wjets.mctperp, wjets.weight, "wjets_template_"+ch, bins, histrange, True)
        # systematic on w+jets template
        rhist = R.TH1D("wjets_syst_"+ch, "wjets_syst_"+ch, bins, histrange[0], histrange[1])
        for i in xrange(bins):
            if templates['wjets_'+ch].GetBinContent(i+1) > 0: #only do non-zero bins
                rhist.SetBinContent(i+1, 0.5) # 50% systematic
        templates['wjets_syst_'+ch] = rhist

        vv = mcvv[selvv['sig_'+ch]]
        templates['vv_'+ch] = rootutils.create_TH1(vv.mctperp, vv.weight, "vv_template_"+ch, bins, histrange, True)

        z = mcz[selz['sig_'+ch]]
        templates['z_'+ch] = rootutils.create_TH1(z.mctperp, z.weight, "z_template_"+ch, bins, histrange, True)

        if ch == 'sf':
            # systematic on Z monte carlo
            mc_onz = mc[smc['z_ctrl_sf']]
            data_onz = data[sd['z_ctrl_sf']]

            mc_hist, mc_edges = np.histogram(mc_onz.mctperp, weights=mc_onz.weight, bins=bins, range=histrange, normed=True)
            d_hist, d_edges = np.histogram(data_onz.mctperp, weights=data_onz.weight, bins=bins, range=histrange, normed=True)

            err = abs(mc_hist[:11]-d_hist[:11])

            # make a TH1 out of it
            rhist = R.TH1D("z_syst", "z_syst", bins, histrange[0], histrange[1])
            for i, val in enumerate(err):
                rhist.SetBinContent(i+1, val)

            templates['z_syst'] = rhist

        n_1tag = sum(mc[smc['1tag_ctrl_'+ch] & (mc.mctperp>5.)].weight)
        n_2tag = sum(mc[smc['2tag_ctrl_'+ch] & (mc.mctperp>5.)].weight)
        eff = 2.*n_2tag/(n_1tag+2*n_2tag)*1.06
        ntop_pred = n_1tag/2.*(1-eff)/eff

        constraints[ch] = R.TVectorD(1)
        constraints[ch][0] = ntop_pred

    for k in templates.keys():
        templates[k].Write()

    for k in constraints.keys():
        constraints[k].Write("ntop_"+k)

    rfile.Close()

def create_data_file(filename="data.root", bins=19, histrange=(10,200)):

    rfile = R.TFile(filename, "RECREATE")
    rfile.cd()

    for ch in channels:
        d = data[sd['sig_'+ch]]
        template = rootutils.create_TH1(d.mctperp, d.weight, "data_"+ch, bins, histrange)
        template.Write()

    rfile.Close()

def create_signal_file(input_file, out_filename, hist_filename, xsec_filename, xsec_multiplier=1., bins=19, histrange=(10,200) ):

    sms_file = HDFStore(input_file)
    sms = sms_file['data']
    sel_sms = selection.get_samples(sms)
    mumu_high_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) > 1.)
    mumu_low_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) < 1.)

    sms.weight *= (sel_sms['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_sms.astype(float)*mumu_high_eta_trigger_eff
                  +mumu_low_eta_sms.astype(float)*mumu_low_eta_trigger_eff + sel_sms['opposite_sign_emu'].astype(float)*emu_trigger_eff)

    # check to see if the file exists, since ROOT will happily continue along with a non-existent file
    if not os.path.exists(hist_filename):
        raise IOError(hist_filename+" does not exist.")
    nevents_file = R.TFile(hist_filename)
    nevents_hist = nevents_file.Get("hist/NEvents_histo")

    xsec_dict = load_xsec(xsec_filename)

    rfile = R.TFile(out_filename, "RECREATE")

    templates = {}
    jes_up_templates = {}
    jes_down_templates = {}
    # create a different template for each mass point
    # templates are scaled to number of exected events given the reference cross-section
    groups = sms.groupby(['mass1', 'mass2'])
    for name, data in groups:
        m1, m2 = name
        sel = selection.get_samples(data)

        events_per_point = nevents_hist.GetBinContent(nevents_hist.FindBin(m1, m2))
        try:
            xsec, xsec_err = xsec_dict[float(m1)]
        except KeyError:
            continue

        xsec *= xsec_multiplier

        for ch in channels:
            mass_point = data[sel['sig_'+ch]]
            mass_point_jes_up = data[sel['sig_scaleup_'+ch]]
            mass_point_jes_down = data[sel['sig_scaledown_'+ch]]

            if mass_point[(mass_point.mctperp > histrange[0]) & (mass_point.mctperp < histrange[1])].mctperp.count() == 0: continue

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

    create_data_file("data.root", bins, histrange)
    create_template_file("templates.root", bins, histrange)
    create_signal_file(args['<signal_file>'], args['<signal_output>'], args['<nevents_file>'], args['<xsec_file>'],
                       float(args['--xsec_multiplier']), bins, histrange)

