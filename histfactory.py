#! /usr/bin/env python

"""Create the model

Usage:
    histfactory.py <signal_file> <template_file> <mass_file> [-h] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    --channels=<c1,c2> Channels to use [default: of,sf]

"""
import ROOT as R
from collections import defaultdict
backgrounds = ['top', 'vv', 'wjets', 'z']


def create_histfactory(template_file, signal_file, m1, m2, channels, data_file_name="data.root", data_prefix="data"):
    prefix = "limits/"+signal_file[:-5]+"_{0}_{1}".format(m1, m2)

    meas = R.RooStats.HistFactory.Measurement("meas", "meas")

    meas.SetOutputFilePrefix(prefix)
    meas.SetPOI("sig_strength")
    # meas.AddConstantParam("n_of_top")
    # meas.AddConstantParam("n_sf_top")

    temp_file = R.TFile("templates.root")

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
        signal = R.RooStats.HistFactory.Sample("signal_"+ch, "sms_template_{0}_{1}_{2}".format(ch, m1, m2), signal_file)
        signal.SetNormalizeByTheory(True)
        signal.AddNormFactor("sig_strength", 1., 0., 1000.)
        signal.ActivateStatError()
        signal.AddOverallSys("trigger", 0.95, 1.05)
        signal.AddOverallSys("id_and_selection", 0.98, 1.02)
        signal.AddOverallSys("b_veto", 0.94, 1.06)
        signal.AddHistoSys("jes", "sms_template_jes_down_{0}_{1}_{2}".format(ch, m1, m2), signal_file, "",
                           "sms_template_jes_up_{0}_{1}_{2}".format(ch, m1, m2), signal_file, "")
        channel_confs[ch].AddSample(signal)

        # add the background samples
        for bkg in backgrounds:
            template = R.RooStats.HistFactory.Sample("{0}_{1}".format(bkg, ch), "{0}_template_{1}".format(bkg, ch), "templates.root")
            template.SetNormalizeByTheory(False)
            template.ActivateStatError()
            template.AddNormFactor("n_{0}_{1}".format(ch, bkg), 2000, 0, 10000)

            if bkg == 'vv':
                template.AddHistoSys('WW_norm', "vv_syst_WW_"+ch+"Up", "templates.root", "", "vv_syst_WW_"+ch+"Down", "templates.root", "")
                template.AddHistoSys('WZ_norm', "vv_syst_WZ_"+ch+"Up", "templates.root", "", "vv_syst_WZ_"+ch+"Down", "templates.root", "")
                template.AddHistoSys('ZZ_norm', "vv_syst_ZZ_"+ch+"Up", "templates.root", "", "vv_syst_ZZ_"+ch+"Down", "templates.root", "")
                template.AddHistoSys('VVV_norm', "vv_syst_VVV_"+ch+"Up", "templates.root", "", "vv_syst_VVV_"+ch+"Down", "templates.root", "")

            if bkg == 'z':
                template.AddShapeSys("z_syst_"+ch, 0, "z_syst", "templates.root")
            if bkg == 'wjets':
                template.AddShapeSys("wjets_syst_"+ch, 0, "wjets_syst_"+ch, "templates.root")


            samples[ch][bkg] = template
            channel_confs[ch].AddSample(samples[ch][bkg])

        meas.AddChannel(channel_confs[ch])

    meas.CollectHistograms()

    R.RooStats.HistFactory.MakeModelAndMeasurementFast(meas)

    print "Saved to", prefix

if __name__ == '__main__':
    from docopt import docopt
    import json

    args = docopt(__doc__)

    sig_file = args['<signal_file>']


    with open(args['<mass_file>']) as f:
        masses = json.load(f)

    try:
        chans = args['--channels'].split(",")
    except AttributeError:
        chans = ['of', 'sf']

    for m1,m2 in masses:
        # try:
        create_histfactory(args['<template_file>'], args['<signal_file>'], int(m1), int(m2), chans)
        break
        # except:
            # continue


