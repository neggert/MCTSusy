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
                rhist.SetBinContent(i+1, 0.3) # 50% systematic
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

            err = abs(mc_hist[:11]-d_hist[:11])/d_hist[:11]

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

def create_combined_file(sms_file, out_filename, nevents_filename, xsec_filename, xsec_multiplier=1., bins=19, histrange=(10,200) ):

    sms_file = HDFStore(sms_file)
    sms = sms_file['data']
    sel_sms = selection.get_samples(sms)
    mumu_high_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) > 1.)
    mumu_low_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) < 1.)

    weight = (sel_sms['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_sms.astype(float)*mumu_high_eta_trigger_eff
                  +mumu_low_eta_sms.astype(float)*mumu_low_eta_trigger_eff + sel_sms['opposite_sign_emu'].astype(float)*emu_trigger_eff)
    weight.name="weight"
    sms = sms.join(weight)

    mcvv = mc[(mc.mc_cat=='WW') | (mc.mc_cat=='ZZ') | (mc.mc_cat=='WZ')]
    mcz = mc[mc.mc_cat=='DY']
    selvv = selection.get_samples( mcvv, 100.)
    selz = selection.get_samples( mcz, 100.)

    # check to see if the file exists, since ROOT will happily continue along with a non-existent file
    if not os.path.exists(nevents_filename):
        raise IOError(nevents_filename+" does not exist.")
    nevents_file = R.TFile(nevents_filename)
    nevents_hist = nevents_file.Get("ScanValidator/NEvents_histo")

    xsec_dict = load_xsec(xsec_filename)

    # rfile = R.TFile(out_filename, "RECREATE")
    # rfile.cd()

    # create a different template for each mass point
    # templates are scaled to number of exected events given the reference cross-section
    groups = sms.groupby(['mass1', 'mass2'])
    for name, sms_data in groups:
        dir_name = "_".join(map(str, map(int, name)))
        ws = R.RooWorkspace(dir_name, dir_name)
        ws.factory("x[{0},{1}]".format(*histrange));
        obs = R.RooArgList(ws.var("x"))
        obs_set = R.RooArgSet(ws.var("x"))
        m1, m2 = name
        sel = selection.get_samples(sms_data)

        events_per_point = nevents_hist.GetBinContent(nevents_hist.FindBin(m1, m2))
        try:
            xsec, xsec_err = xsec_dict[float(m1)]
        except KeyError:
            continue

        xsec *= xsec_multiplier

        for ch in channels:
            # signal templates
            mass_point = sms_data[sel['sig_'+ch]]
            mass_point_jes_up = sms_data[sel['sig_scaleup_'+ch]]
            mass_point_jes_down = sms_data[sel['sig_scaledown_'+ch]]

            if mass_point[(mass_point.mctperp > histrange[0]) & (mass_point.mctperp < histrange[1])].mctperp.count() == 0: continue

            h = rootutils.create_TH1(mass_point.mctperp, mass_point.weight,
                                                                         "sms_template",
                                                                          bins, histrange, False)
            hup = rootutils.create_TH1(mass_point_jes_up.mctperp_up, mass_point_jes_up.weight,
                                                                         "sms_template_jesUp",
                                                                          bins, histrange, False)
            hdown = rootutils.create_TH1(mass_point_jes_down.mctperp_down, mass_point_jes_down.weight,
                                                                         "sms_template_jesDown",
                                                                          bins, histrange, False)
            h.Scale(xsec*lumi/events_per_point)
            hup.Scale(xsec*lumi/events_per_point)
            hdown.Scale(xsec*lumi/events_per_point)
            h_rdh = R.RooDataHist(ch+"_sms_template",ch+"_sms_template", obs, h, 1.0)
            getattr(ws, "import")(h_rdh)
            hup_rdh = R.RooDataHist(ch+"_sms_template_jesUp",ch+"_sms_template_jesUp", obs, hup)
            getattr(ws, "import")(hup_rdh)
            hdown_rdh = R.RooDataHist(ch+"_sms_template_jesDown",ch+"_sms_template_jesDown", obs, hdown)
            getattr(ws, "import")(hdown_rdh)


            # background templates
            top = data[sd['top_ctrl_'+ch]]
            top_template  = rootutils.create_TH1(top.mctperp, top.weight, "top_template", bins, histrange, True)
            top_rdh = R.RooDataHist(ch+"_top_template",ch+"_top_template", obs, top_template)
            getattr(ws, "import")(top_rdh)
            top_rhf = R.RooHistPdf(ch+"_top_template_func", ch+"_top_template_func", obs_set, ws.obj(ch+"_top_template"))
            getattr(ws, "import")(top_rhf)
            ws.factory(ch+"_top_template_pdf_norm[0,10000]")
            top_pdf = R.RooExtendPdf(ch+"_top_template_pdf", ch+"_top_template_pdf", ws.obj(ch+"_top_template_func"),
                                     ws.var(ch+"_top_template_pdf_norm"))
            getattr(ws, "import")(top_pdf)

            wjets = data[sd['wjets_ctrl_'+ch]]
            wjets_template = rootutils.create_TH1(wjets.mctperp, wjets.weight, "wjets_template", bins, histrange, True)
            wjets_rdh = R.RooDataHist(ch+"_wjets_template",ch+"_wjets_template", obs, wjets_template)
            getattr(ws, "import")(wjets_rdh)
            wjets_rhf = R.RooHistPdf(ch+"_wjets_template_func", ch+"_wjets_template_func", obs_set, ws.obj(ch+"_wjets_template"))
            getattr(ws, "import")(wjets_rhf)
            ws.factory(ch+"_wjets_template_pdf_norm[0,10000]")
            wjets_pdf = R.RooExtendPdf(ch+"_wjets_template_pdf", ch+"_wjets_template_pdf", ws.obj(ch+"_wjets_template_func"),
                                     ws.var(ch+"_wjets_template_pdf_norm"))
            getattr(ws, "import")(wjets_pdf)
            # systematic on w+jets template
            rhist = R.TH1D("wjets_syst_"+ch, "wjets_syst", bins, histrange[0], histrange[1])
            for i in xrange(bins):
                if wjets_template.GetBinContent(i+1) > 0: #only do non-zero bins
                    rhist.SetBinContent(i+1, 0.3) # 50% systematic
            rhist_rdh = R.RooDataHist(ch+"_wjets_syst",ch+"_wjets_syst", obs, rhist)
            getattr(ws,"import")(rhist_rdh)

            vv = mcvv[selvv['sig_'+ch]]
            vv_template = rootutils.create_TH1(vv.mctperp, vv.weight, "vv_template", bins, histrange, True)
            vv_rdh = R.RooDataHist(ch+"_vv_template",ch+"_vv_template", obs, vv_template)
            getattr(ws, "import")(vv_rdh)
            vv_rhf = R.RooHistPdf(ch+"_vv_template_func", ch+"_vv_template_func", obs_set, ws.obj(ch+"_vv_template"))
            getattr(ws, "import")(vv_rhf)
            ws.factory(ch+"_vv_template_pdf_norm[0,10000]")
            vv_pdf = R.RooExtendPdf(ch+"_vv_template_pdf", ch+"_vv_template_pdf", ws.obj(ch+"_vv_template_func"),
                                     ws.var(ch+"_vv_template_pdf_norm"))
            getattr(ws, "import")(vv_pdf)

            z = mcz[selz['sig_'+ch]]
            z_template = rootutils.create_TH1(z.mctperp, z.weight, "z_template", bins, histrange, True)
            z_rdh = R.RooDataHist(ch+"_z_template",ch+"_z_template", obs, z_template)
            getattr(ws, "import")(z_rdh)
            z_rhf = R.RooHistPdf(ch+"_z_template_func", ch+"_z_template_func", obs_set, ws.obj(ch+"_z_template"))
            getattr(ws, "import")(z_rhf)
            ws.factory(ch+"_z_template_pdf_norm[0,10000]")
            z_pdf = R.RooExtendPdf(ch+"_z_template_pdf", ch+"_z_template_pdf", ws.obj(ch+"_z_template_func"),
                                     ws.var(ch+"_z_template_pdf_norm"))
            getattr(ws, "import")(z_pdf)

            if ch == 'sf':
                # systematic on Z monte carlo
                mc_onz = mc[smc['z_ctrl_sf']]
                data_onz = data[sd['z_ctrl_sf']]

                mc_hist, mc_edges = np.histogram(mc_onz.mctperp, weights=mc_onz.weight, bins=bins, range=histrange, normed=True)
                d_hist, d_edges = np.histogram(data_onz.mctperp, weights=data_onz.weight, bins=bins, range=histrange, normed=True)

                err = abs(mc_hist[:11]-d_hist[:11])/d_hist[:11]

                # make a TH1 out of it
                rhist = R.TH1D("z_syst", "z_syst", bins, histrange[0], histrange[1])
                for i, val in enumerate(err):
                    rhist.SetBinContent(i+1, val)

                rhist_rdh = R.RooDataHist(ch+"_z_syst",ch+"_z_syst", obs, rhist)
                getattr(ws,"import")(rhist_rdh)
            # data
            d = data[sd['sig_'+ch]]
            template = rootutils.create_TH1(d.mctperp, d.weight, "data_obs", bins, histrange)
            data_rdh = R.RooDataHist(ch+"_data_obs",ch+"_data_obs", obs, template)
            getattr(ws, "import")(data_rdh)

        ws.writeToFile(out_filename, False)
        break

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)

    print(args)

    bins = int(args['--bins'])
    histrange = (float(args['--low']), float(args['--high']))

    create_combined_file(args['<signal_file>'], args['<signal_output>'], args['<nevents_file>'], args['<xsec_file>'],
                       float(args['--xsec_multiplier']), bins, histrange)

