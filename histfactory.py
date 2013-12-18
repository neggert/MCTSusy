#! /usr/bin/env python

"""Create the model

Usage:
    histfactory.py <template_file> [<signal_file>] [<output>] [-h] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    --channels=<c1,c2> Channels to use [default: of,sf]

"""
import ROOT as R
from collections import defaultdict
import sys
import re
backgrounds = ['top', 'vv', 'wjets', 'z']


def create_histfactory(template_file, channels, data_file_name="data.root", signal_file=None, m1=0, m2=0):
    data_prefix="data"
    meas = R.RooStats.HistFactory.Measurement("meas", "meas")

    if signal_file:
        prefix = "limits/"+signal_file[:-5]+"_{0}_{1}".format(m1, m2)
        meas.SetOutputFilePrefix(prefix)
    else:
        meas.SetOutputFilePrefix("limits/chi_bkg_only")
    # meas.AddConstantParam("n_of_top")
    # meas.AddConstantParam("n_sf_top")
    meas.SetPOI("sig_strength") 

    temp_file = R.TFile(template_file)

    meas.SetLumi(1.0)
    meas.SetLumiRelErr(0.04)
    meas.SetExportOnly(True)

    channel_confs = {}
    samples = defaultdict(dict)

    for ch in channels:
        channel_confs[ch] = R.RooStats.HistFactory.Channel(ch)
        channel_confs[ch].SetData(data_prefix+"_"+ch, data_file_name)
        channel_confs[ch].SetStatErrorConfig(0.1, "Poisson")

        # signal sample
        if signal_file:
            signal = R.RooStats.HistFactory.Sample("signal_"+ch, "sms_template_{0}_{1}_{2}".format(ch, m1, m2), signal_file)
            signal.SetNormalizeByTheory(True)
            signal.AddNormFactor("sig_strength", 1., 0., 1000.)
            signal.ActivateStatError()
            signal.AddOverallSys("trigger", 0.95, 1.05)
            signal.AddOverallSys("id_and_selection", 0.976, 1.024)
            signal.AddOverallSys("b_veto", 0.95, 1.05)
            signal.AddHistoSys("jes", "sms_template_jes_down_{0}_{1}_{2}".format(ch, m1, m2), signal_file, "",
                               "sms_template_jes_up_{0}_{1}_{2}".format(ch, m1, m2), signal_file, "")
            channel_confs[ch].AddSample(signal)

        # add the background samples
        for bkg in backgrounds:
            template = R.RooStats.HistFactory.Sample("{0}_{1}".format(bkg, ch), "{0}_template_{1}".format(bkg, ch), "templates.root")
            template.SetNormalizeByTheory(False)
            template.ActivateStatError()
            if bkg=="top":
                template.AddNormFactor("n_top_sf".format(ch, bkg), 1, 0, 10000)
                if ch=="of":
                    top_ratio_val = temp_file.Get("top_ratio")[0]
                    template.AddNormFactor("n_top_of_scale", top_ratio_val, top_ratio_val, top_ratio_val, True)
                    template.AddOverallSys("top_ratio", .9, 1.1)
                else:
                    template.AddOverallSys("top_ratio", 1.1, 0.9)
            elif bkg=="vv":
                template.AddNormFactor("n_vv_sf".format(ch, bkg), 1, 0, 10000)
                if ch=="of":
                    vv_ratio_val = temp_file.Get("vv_ratio")[0]
                    template.AddNormFactor("n_vv_of_scale", vv_ratio_val, vv_ratio_val, vv_ratio_val, True)
                    template.AddOverallSys("vv_ratio", .9, 1.1)
                else:
                    template.AddOverallSys("vv_ratio", 1.1, 0.9)
            else:
                template.AddNormFactor("n_{0}_{1}".format(ch, bkg), 1, 0, 10000)

            if bkg == 'vv':
                template.AddShapeSys("ww_syst_"+ch, 0, "ww_syst_"+ch, "templates.root")
                template.AddHistoSys('WW_norm', "vv_syst_WW_"+ch+"Up", "templates.root", "", "vv_syst_WW_"+ch+"Down", "templates.root", "")
                template.AddHistoSys('WZ_norm', "vv_syst_WZ_"+ch+"Up", "templates.root", "", "vv_syst_WZ_"+ch+"Down", "templates.root", "")
                template.AddHistoSys('ZZ_norm', "vv_syst_ZZ_"+ch+"Up", "templates.root", "", "vv_syst_ZZ_"+ch+"Down", "templates.root", "")
                template.AddHistoSys('VVV_norm', "vv_syst_VVV_"+ch+"Up", "templates.root", "", "vv_syst_VVV_"+ch+"Down", "templates.root", "")
                template.AddHistoSys('HWW_norm', "vv_syst_HWW_"+ch+"Up", "templates.root", "", "vv_syst_HWW_"+ch+"Down", "templates.root", "")


            if bkg == 'z':
                template.AddShapeSys("z_syst_"+ch, 0, "z_syst", "templates.root")
            if bkg == 'wjets':
                template.AddShapeSys("wjets_syst_"+ch, 0, "wjets_syst_"+ch, "templates.root")


            samples[ch][bkg] = template
            channel_confs[ch].AddSample(samples[ch][bkg])

        meas.AddChannel(channel_confs[ch])

    meas.CollectHistograms()

    ws = R.RooStats.HistFactory.MakeModelAndMeasurementFast(meas)

    # top_ratio_val = temp_file.Get("top_ratio")[0]
    # ws.factory('expr::top_ratio("n_of_top/n_sf_top", n_of_top, n_sf_top)')
    # ws.factory('RooGaussian::top_ratio_constraint(top_ratio, nom_top_ratio[{0}], {1})'.format(top_ratio_val, top_ratio_val*0.1))
    # vv_ratio_val = temp_file.Get("vv_ratio")[0]
    # ws.factory('expr::vv_ratio("n_of_vv/n_sf_vv", n_of_vv, n_sf_vv)')
    # ws.factory('RooGaussian::vv_ratio_constraint(vv_ratio, nom_vv_ratio[{0}], {1})'.format(vv_ratio_val, vv_ratio_val*0.1))
    # ws.factory('PROD:constrPdf(simPdf, top_ratio_constraint, vv_ratio_constraint)')

    # model = ws.obj("ModelConfig")
    # model.SetPdf(ws.obj("constrPdf"))

    # # ws.Print()

    # rfile = R.TFile(prefix+"_constrained.root", "RECREATE")
    # ws.Write()
    # rfile.Close()

if __name__ == '__main__':
    from docopt import docopt
    import json

    args = docopt(__doc__)
    sys.argv = [sys.argv[0], "-b"]

    sig_file = args['<signal_file>']

    try:
        chans = args['--channels'].split(",")
    except AttributeError:
        chans = ['of', 'sf']

    if sig_file:
        m1, m2 = re.search("_(\d*?)_(\d*?)_", args['<output>']).groups()[:2]

        create_histfactory(args['<template_file>'], chans, "data.root", args['<signal_file>'], int(m1), int(m2))

    else:
        create_histfactory(args['<template_file>'], chans, "data.root")



