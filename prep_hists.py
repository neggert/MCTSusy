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

import sys
sys.path.append("import")
import CMSPyLibs.general_calc as general_calc

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

    mcvv = mc[(mc.mc_cat=='WW') | (mc.mc_cat=='ZZ') | (mc.mc_cat=='WZ') | (mc.mc_cat=='VVV') | (mc.mc_cat=='HWW')]
    mcz = mc[mc.mc_cat=='DY']
    selvv = selection.get_samples( mcvv, 100.)
    selz = selection.get_samples( mcz, 100.)

    rfile = R.TFile(filename, "RECREATE")
    rfile.cd()

    templates = {}

    constraints = {}

    # prep WW systematics
    three_m1_mc = mc[mc.ThirdLepton & ~np.isnan(mc.pt3) & (mc.metPt > 30)]
    three_m1_mc = three_m1_mc.apply(shuffle_leptons, axis=1)

    three_m1_mc.mctperp = three_m1_mc.apply(recalc_MCT, axis=1)
    three_m1_mc.mctperp *= 80.4/91.2
    three_m1_mc.metPt = three_m1_mc.apply(recalc_MET, axis=1)
    three_m1_mc.ThirdLepton = False
    sthree_m1_mc = selection.get_samples(three_m1_mc)
    three_m1_mc = three_m1_mc[sthree_m1_mc['sig_of']]

    three_m1_data = data[data.ThirdLepton]
    three_m1_data = three_m1_data.apply(shuffle_leptons, axis=1)
    three_m1_data.mctperp = three_m1_data.apply(recalc_MCT, axis=1)
    three_m1_data.mctperp *= 80.4/91.2
    three_m1_data.metPt = three_m1_data.apply(recalc_MET, axis=1)
    three_m1_data.ThirdLepton = False
    sthree_m1_data = selection.get_samples(three_m1_data)
    three_m1_data = three_m1_data[sthree_m1_data['sig_of']]

    ww = mc[smc['sig_of'] & (mc.mc_cat=="WW")]

    # re-weighting factors from comparing 3-1 lepton region to WW in MC
    ww_hist, ww_bins = np.histogram(ww.mctperp, weights=ww.weight, bins=14, range=(10,290), normed=True)
    three_m1_mc_hist, _ = np.histogram(three_m1_mc.mctperp, weights=three_m1_mc.weight, bins=14, range=(10,290), normed=True)
    three_m1_reweighting = ww_hist/three_m1_mc_hist

    def get_weight(mctperp):
        i = np.argmax(np.where(ww_bins<mctperp, ww_bins, 0))
        return three_m1_reweighting[i]

    # re-weight data 3-1 region
    reweighting = three_m1_data.mctperp.apply(get_weight)
    three_m1_data.weight *= reweighting

    # derive systematic from comparison of 3-1 data to MC
    # systematic is the larger of fractional uncertainty on data and discrepancy between data and MC
    three_m1_data_hist, _ = np.histogram(three_m1_data.mctperp, bins=14, range=(10,290)) # unweighted
    three_m1_stat = 1./np.sqrt(three_m1_data_hist) # fractional
    three_m1_data_hist_normed, _ = np.histogram(three_m1_data.mctperp, weights=three_m1_data.weight, bins=14, range=(10,290), normed=True)
    ww_discrepancy = abs(three_m1_data_hist_normed-ww_hist)/ww_hist

    ww_systematic = np.max((three_m1_stat, ww_discrepancy), axis=0)
    ww_systematic[np.isinf(ww_systematic)] = 1.

    for ch in channels:
        top = data[sd['top_ctrl_'+ch]]
        templates['top_'+ch] = rootutils.create_TH1(top.mctperp, top.weight, "top_template_"+ch, bins, histrange, True)

        wjets = data[sd['wjets_ctrl_'+ch]]
        templates['wjets_'+ch] = rootutils.create_TH1(wjets.mctperp, wjets.weight, "wjets_template_"+ch, bins, histrange, True)
        # systematic on w+jets template
        rhist = R.TH1D("wjets_syst_"+ch, "wjets_syst_"+ch, bins, histrange[0], histrange[1])
        for i in xrange(bins):
            # if templates['wjets_'+ch].GetBinContent(i+1) > 0: #only do non-zero bins
            rhist.SetBinContent(i+1, 0.3) # 50% systematic
        templates['wjets_syst_'+ch] = rhist

        vv = mcvv[selvv['sig_'+ch]]
        templates['vv_'+ch] = rootutils.create_TH1(vv.mctperp, vv.weight, "vv_template_"+ch, bins, histrange, True)

        ww_hist, template_bins = np.histogram(ww.mctperp, weights=ww.weight, bins=bins, range=histrange)
        vv_hist, template_bins = np.histogram(vv.mctperp, weights=vv.weight, bins=bins, range=histrange)
        ww_frac = ww_hist/vv_hist
        ww_frac[np.isnan(ww_frac)] = 0.
        bin_centers = (template_bins[1:]+template_bins[:-1])/2
        rhist = R.TH1D("ww_syst_"+ch, "ww_syst_"+ch, bins, histrange[0], histrange[1])

        for i, b in enumerate(bin_centers):
            if b > ww_bins[-1]:
                syst = 0.
            else:
                j = np.argmax(np.where(ww_bins<b, ww_bins, 0))
                syst = ww_systematic[j]*ww_frac[i]
                if np.isnan(syst): syst = 0.
            rhist.SetBinContent(i+1, syst)
        templates['ww_syst_'+ch] = rhist

        # shape systematic
        weights = vv.weight*(1+(vv.mc_cat=="WW")*0.10)
        templates['vv_syst_WW_'+ch+'Up'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_WW_"+ch+"Up", bins, histrange, True)
        weights = vv.weight*(1-(vv.mc_cat=="WW")*0.10)
        templates['vv_syst_WW_'+ch+'Down'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_WW_"+ch+"Down", bins, histrange, True)

        weights = vv.weight*(1+(vv.mc_cat=="WZ")*0.10)
        templates['vv_syst_WZ_'+ch+'Up'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_WZ_"+ch+"Up", bins, histrange, True)
        weights = vv.weight*(1-(vv.mc_cat=="WZ")*0.10)
        templates['vv_syst_WZ_'+ch+'Down'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_WZ_"+ch+"Down", bins, histrange, True)

        weights = vv.weight*(1+(vv.mc_cat=="ZZ")*0.10)
        templates['vv_syst_ZZ_'+ch+'Up'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_ZZ_"+ch+"Up", bins, histrange, True)
        weights = vv.weight*(1-(vv.mc_cat=="ZZ")*0.10)
        templates['vv_syst_ZZ_'+ch+'Down'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_ZZ_"+ch+"Down", bins, histrange, True)

        weights = vv.weight*(1+(vv.mc_cat=="VVV")*0.50)
        templates['vv_syst_VVV_'+ch+'Up'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_VVV_"+ch+"Up", bins, histrange, True)
        weights = vv.weight*(1-(vv.mc_cat=="VVV")*0.50)
        templates['vv_syst_VVV_'+ch+'Down'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_VVV_"+ch+"Down", bins, histrange, True)

        weights = vv.weight*(1+(vv.mc_cat=="HWW")*0.50)
        templates['vv_syst_HWW_'+ch+'Up'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_HWW_"+ch+"Up", bins, histrange, True)
        weights = vv.weight*(1-(vv.mc_cat=="HWW")*0.50)
        templates['vv_syst_HWW_'+ch+'Down'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_HWW_"+ch+"Down", bins, histrange, True)

        z = mcz[selz['sig_'+ch]]
        templates['z_'+ch] = rootutils.create_TH1(z.mctperp, z.weight, "z_template_"+ch, bins, histrange, True)

        if ch == 'sf':
            # systematic on Z monte carlo
            mc_onz = mc[smc['z_ctrl_sf']]
            data_onz = data[sd['z_ctrl_sf']]

            mc_hist, mc_edges = np.histogram(mc_onz.mctperp, weights=mc_onz.weight, bins=bins, range=histrange, normed=True)
            d_hist, d_edges = np.histogram(data_onz.mctperp, weights=data_onz.weight, bins=bins, range=histrange, normed=True)

            err = np.zeros(bins)
            err[:] = np.max(abs(mc_hist[:10]-d_hist[:10])/d_hist[:10])
            # err[10:] = 0.5

            # make a TH1 out of it
            rhist = R.TH1D("z_syst", "z_syst", bins, histrange[0], histrange[1])
            for i, val in enumerate(err):
                rhist.SetBinContent(i+1, val)

            templates['z_syst'] = rhist

        n_1tag = sum(mc[smc['1tag_ctrl_'+ch] & (mc.mctperp>5.)].weight)
        n_2tag = sum(mc[smc['2tag_ctrl_'+ch] & (mc.mctperp>5.)].weight)
        eff = 2.*n_2tag/(n_1tag+2*n_2tag)*1.06
        ntop_pred = n_1tag/2.*(1-eff)/eff


    constraints['top_ratio'] = R.TVectorD(1)
    constraints['top_ratio'][0] = mc[smc['sig_mct_low_of'] & (mc.mc_cat=="top")].weight.sum()/mc[smc['sig_mct_low_sf'] & (mc.mc_cat=="top")].weight.sum()
    constraints['vv_ratio'] = R.TVectorD(1)
    constraints['vv_ratio'][0] = (mc[smc['sig_mct_low_of'] & (mc.mc_cat=="WW")].weight.sum()/mc[smc['sig_mct_low_sf'] & (mc.mc_cat=="WW")].weight.sum())\
                                 * mcvv[selvv['sig_mct_low_sf']].weight.sum()/mcvv[selvv['sig_mct_low_of']].weight.sum()

    for k in templates.keys():
        templates[k].Write()

    for k in constraints.keys():
        constraints[k].Write(k)

    rfile.Close()

def create_data_file(filename="data.root", bins=19, histrange=(10,200), bootstrap_num=0):

    rfile = R.TFile(filename, "RECREATE")
    rfile.cd()

    bs_label = ""
    for ch in channels:
        if bootstrap_num == 0:
            d = data[sd['sig_'+ch]]
        else:
            rand = np.random.mtrand.RandomState(bootstrap_num*100000+333)
            d = data[sd['sig_'+ch]]
            data_size = len(d.index)
            d = d.ix[d.index[rand.choice(data_size, size=data_size, replace=True)]]
            bs_label = "bs_{0}_".format(bootstrap_num)
        template = rootutils.create_TH1(d.mctperp, d.weight, "data_"+bs_label+ch, bins, histrange)
        template.Write()

    rfile.Close()

def create_data_file_with_signal(sms_filename, xsec_filename, hist_filename, filename="data.root", bins=19, histrange=(10,200),
                                 mass1=400, mass2=100, xsec_multiplier=1.):



    sms_file = HDFStore(sms_filename)
    sms = sms_file['data']
    sel_sms = selection.get_samples(sms)
    mumu_high_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) > 1.)
    mumu_low_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) < 1.)

    weight = (sel_sms['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_sms.astype(float)*mumu_high_eta_trigger_eff
                  +mumu_low_eta_sms.astype(float)*mumu_low_eta_trigger_eff + sel_sms['opposite_sign_emu'].astype(float)*emu_trigger_eff)
    weight.name="weight"
    sms = sms.join(weight)

    sms_data = sms[(sms.mass1==mass1) & (sms.mass2==mass2)]
    sel = selection.get_samples(sms_data)


    # check to see if the file exists, since ROOT will happily continue along with a non-existent file
    if not os.path.exists(hist_filename):
        raise IOError(hist_filename+" does not exist.")
    nevents_file = R.TFile(hist_filename)
    nevents_hist = nevents_file.Get("ScanValidator/NEvents_histo")

    xsec_dict = load_xsec(xsec_filename)

    events_per_point = nevents_hist.GetBinContent(nevents_hist.FindBin(mass1, mass2))
    xsec, xsec_err = xsec_dict[float(mass1)]

    xsec *= xsec_multiplier

    rfile = R.TFile(filename, "RECREATE")
    rfile.cd()

    for ch in channels:
        d = data[sd['sig_'+ch]]
        s = sms_data[sel['sig_'+ch]]
        template = rootutils.create_TH1(d.mctperp, d.weight, "data_"+ch, bins, histrange)
        sms_template = rootutils.create_TH1(s.mctperp, s.weight, "sms_"+ch, bins, histrange)
        sms_template.Scale(xsec*lumi/events_per_point)
        for i in xrange(int(sms_template.Integral())):
            template.Fill(sms_template.GetRandom())
        template.Write()

    rfile.Close()

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

    # check to see if the file exists, since ROOT will happily continue along with a non-existent file
    if not os.path.exists(hist_filename):
        raise IOError(hist_filename+" does not exist.")
    nevents_file = R.TFile(hist_filename)
    nevents_hist = nevents_file.Get("ScanValidator/NEvents_histo")

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


calc = general_calc.VarCalculator()

def recalc_MCT(event):
    p4jet1 = R.TLorentzVector(0., 0., 0., 0.)
    p4jet2 = R.TLorentzVector(0., 0., 0., 0.)
    p4l1 = R.TLorentzVector(0., 0., 0., 0.)
    p4l2 = R.TLorentzVector(0., 0., 0., 0.)
    p4l1.SetPtEtaPhiM(event['pt1'], event['eta1'], event['phi1'], 0.)
    p4l2.SetPtEtaPhiM(event['pt2'], event['eta2'], event['phi2'], 0.)

    # Add third lepton into MET
    # note that we don't care about eta or M for these vectors, as they're not used.
    p4met = R.TLorentzVector(0., 0., 0., 0.)
    p4met.SetPtEtaPhiM(event['metPt'], 0., event['metPhi'], 0.)
    p4l3 = R.TLorentzVector(0., 0., 0., 0.)
    p4l3.SetPtEtaPhiM(event['pt3'], 0., event['phi3'], 0.)

    p4met += p4l3

    calc.setP4s(p4jet1, p4jet2, p4l1, p4l2, p4met)
    return calc.mctPerp_210()

def recalc_MET(event):

    # Add third lepton into MET
    # note that we don't care about eta or M for these vectors, as they're not used.
    p4met = R.TLorentzVector(0., 0., 0., 0.)
    p4met.SetPtEtaPhiM(event['metPt'], 0., event['metPhi'], 0.)
    p4l3 = R.TLorentzVector(0., 0., 0., 0.)
    p4l3.SetPtEtaPhiM(event['pt3'], 0., event['phi3'], 0.)

    p4met += p4l3

    return p4met.Pt()

# do some shuffling to recover as many events as possible
def shuffle_leptons(event):
    # make an OF pair if possible
    if abs(event['pdg1']) == abs(event['pdg2']):
        if abs(event['pdg1']) != abs(event['pdg3']):
            event['pdg2'], event['pdg3'] = event['pdg3'], event['pdg2']
            event['pt2'], event['pt3'] = event['pt3'], event['pt2']
            event['eta2'], event['eta3'] = event['eta3'], event['eta2']
            event['phi2'], event['phi3'] = event['phi3'], event['phi2']

        elif abs(event['pdg2']) != abs(event['pdg3']):
            event['pdg1'], event['pdg3'] = event['pdg3'], event['pdg1']
            event['pt1'], event['pt3'] = event['pt3'], event['pt1']
            event['eta1'], event['eta3'] = event['eta3'], event['eta1']
            event['phi1'], event['phi3'] = event['phi3'], event['phi1']

    return event

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

