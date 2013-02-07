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
    backgrounds = ['of', 'vv', 'wjets', 'z']

    m1 = int(m1)
    m2 = int(m2)

    meas = R.RooStats.HistFactory.Measurement("meas", "meas")

    meas.SetOutputFilePrefix(prefix)
    meas.SetPOI("sig_strength")
    # meas.AddConstantParam("n_of_top")
    # meas.AddConstantParam("n_sf_top")

    meas.SetLumi(1.0)
    meas.SetLumiRelErr(0.04)
    meas.SetExportOnly(True)

    channel_conf = None
    samples = {}

    channel_conf = R.RooStats.HistFactory.Channel('sf')
    channel_conf.SetData("data_sf", "data.root")
    channel_conf.SetStatErrorConfig(0.1, "Poisson")

    # signal sample
    signal = R.RooStats.HistFactory.Sample("signal_sf", "sms_template_sf_{0}_{1}".format(m1, m2), signal_file)
    signal.SetNormalizeByTheory(True)
    signal.AddNormFactor("sig_strength", 1., 0., 100.)
    signal.ActivateStatError()
    signal.AddOverallSys("trigger", 0.95, 1.05)
    signal.AddOverallSys("id_and_selection", 0.98, 1.02)
    signal.AddOverallSys("b_veto", 0.94, 1.06)
    signal.AddHistoSys("jes", "sms_template_jes_down_sf_{0}_{1}".format(m1, m2), signal_file, "",
                       "sms_template_jes_up_sf_{0}_{1}".format(m1, m2), signal_file, "")
    channel_conf.AddSample(signal)

    # add the background samples
    for bkg in backgrounds:
        template = R.RooStats.HistFactory.Sample("{0}".format(bkg), "{0}_template".format(bkg), "templates2.root")
        template.SetNormalizeByTheory(False)
        template.ActivateStatError()
        template.AddNormFactor("n_{0}".format(bkg), 2000, 0, 5000)

        if bkg == 'z':
            template.AddShapeSys("z_syst", 0, "z_syst", "templates2.root")
        if bkg == 'wjets':
            template.AddShapeSys("wjets_syst", 0, "wjets_syst", "templates2.root")

        samples[bkg] = template
        channel_conf.AddSample(samples[bkg])

    meas.AddChannel(channel_conf)

    meas.CollectHistograms()

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
        try:
            create_histfactory(args['<template_file>'], args['<signal_file>'], int(m1), int(m2), chans)
        except:
            continue

