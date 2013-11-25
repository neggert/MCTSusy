#! /usr/bin/env python

"""Create the model

Usage:
    histfactory.py <signal_file> <template_file> <mass_file> [-h] 

Options:
    -h --help        Show this screen.

"""
import ROOT as R

backgrounds = ['of', 'vv', 'wjets', 'z']

def create_histfactory(template_file, signal_file, m1, m2, data_file_name="data.root", data_prefix="data"):
    prefix = "limits/"+signal_file[:-5]+"_{0}_{1}".format(m1, m2)

    meas = R.RooStats.HistFactory.Measurement("meas", "meas")

    meas.SetOutputFilePrefix(prefix)
    meas.SetPOI("sig_strength")

    meas.SetLumi(1.0)
    meas.SetLumiRelErr(0.04)
    meas.SetExportOnly(True)

    samples = {}

    channel_conf = R.RooStats.HistFactory.Channel('sf')
    channel_conf.SetData("data_sf", "data.root")
    channel_conf.SetStatErrorConfig(0.01, "Poisson")

    # signal sample
    signal = R.RooStats.HistFactory.Sample("signal_sf", "sms_template_sf_{0}_{1}".format(m1, m2), signal_file)
    signal.SetNormalizeByTheory(True)
    signal.AddNormFactor("sig_strength", 1., 0., 100.)
    signal.ActivateStatError()
    signal.AddOverallSys("trigger", 0.95, 1.05)
    signal.AddOverallSys("id_and_selection", 0.98, 1.02)
    signal.AddOverallSys("b_veto", 0.95, 1.05)
    signal.AddHistoSys("jes", "sms_template_jes_down_sf_{0}_{1}".format(m1, m2), signal_file, "",
                       "sms_template_jes_up_sf_{0}_{1}".format(m1, m2), signal_file, "")
    channel_conf.AddSample(signal)

    # add the background samples
    for bkg in backgrounds:
        template = R.RooStats.HistFactory.Sample("{0}".format(bkg), "{0}_template".format(bkg), template_file)
        template.SetNormalizeByTheory(False)
        template.ActivateStatError()
        template.AddNormFactor("n_{0}".format(bkg), 2000, 0, 10000)

        if bkg == 'z':
            template.AddShapeSys("z_syst", 0, "z_syst", template_file)
        if bkg == 'wjets':
            template.AddShapeSys("wjets_syst", 0, "wjets_syst", template_file)
        if bkg == "vv":
            template.AddHistoSys('WZ_norm', "vv_syst_WZ_Up", "templates2.root", "", "vv_syst_WZ_Down", "templates2.root", "")
            template.AddHistoSys('ZZ_norm', "vv_syst_ZZ_Up", "templates2.root", "", "vv_syst_ZZ_Down", "templates2.root", "")

        samples[bkg] = template
        channel_conf.AddSample(samples[bkg])

    meas.AddChannel(channel_conf)

    meas.CollectHistograms()

    R.RooStats.HistFactory.MakeModelAndMeasurementFast(meas)

if __name__ == '__main__':
    from docopt import docopt
    import json

    args = docopt(__doc__)

    sig_file = args['<signal_file>']


    with open(args['<mass_file>']) as f:
        masses = json.load(f)

    for m1,m2 in masses:
        try:
            create_histfactory(args['<template_file>'], args['<signal_file>'], int(m1), int(m2), "data2.root")
            # break
        except:
            continue

