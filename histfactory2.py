#! /usr/bin/env python

"""Create the model

Usage:
    histfactory2.py <template_file> [<signal_file>] [<mass_file>] [-h] 

Options:
    -h --help        Show this screen.

"""
import ROOT as R

backgrounds = ['of', 'vv', 'wjets', 'z']

def create_histfactory(template_file, signal_file=None, m1=0, m2=0):
    data_prefix="data"
    data_file_name="data.root"
    meas = R.RooStats.HistFactory.Measurement("meas", "meas")

    if signal_file:
        prefix = "limits/"+signal_file[:-5]+"_{0}_{1}".format(m1, m2)
        meas.SetOutputFilePrefix(prefix)
    else:
        meas.SetOutputFilePrefix("limits/slep_bkg_only")
    meas.SetPOI("sig_strength")

    meas.SetLumi(1.0)
    meas.SetLumiRelErr(0.04)
    meas.SetExportOnly(True)

    samples = {}

    channel_conf = R.RooStats.HistFactory.Channel('sf')
    channel_conf.SetData("data_sf", "data.root")
    channel_conf.SetStatErrorConfig(0.01, "Poisson")

    # signal sample
    if signal_file:
        signal = R.RooStats.HistFactory.Sample("signal_sf", "sms_template_sf_{0}_{1}".format(m1, m2), signal_file)
        signal.SetNormalizeByTheory(True)
        signal.AddNormFactor("sig_strength", 1., 0., 100.)
        signal.ActivateStatError()
        signal.AddOverallSys("trigger", 0.95, 1.05)
        signal.AddOverallSys("id_and_selection", 0.976, 1.024)
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
            template.AddHistoSys('WZ_norm', "vv_syst_WZ_Up", template_file, "", "vv_syst_WZ_Down", template_file, "")
            template.AddHistoSys('ZZ_norm', "vv_syst_ZZ_Up", template_file, "", "vv_syst_ZZ_Down", template_file, "")

        samples[bkg] = template
        channel_conf.AddSample(samples[bkg])

    meas.AddChannel(channel_conf)

    meas.CollectHistograms()

    R.RooStats.HistFactory.MakeModelAndMeasurementFast(meas)

if __name__ == '__main__':
    from docopt import docopt
    import json

    args = docopt(__doc__)

    for m1,m2 in masses:
        try:
            create_histfactory(args['<template_file>'], args['<signal_file>'], int(m1), int(m2), "data2.root")
            # break
        except:
            continue

