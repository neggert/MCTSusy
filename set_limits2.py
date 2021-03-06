#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    set_limits2.py <signal_file> <mass1> <mass2> [-hae] [--ncpu=<c>] [--channels=<c1,c2>]
    set_limits2.py batch <signal_file> <mass_file> <jobnum> <output_file> [-aeh] [--ncpu=<c>] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    -a --asymptotic  Run asymptotic limits
    -e --expected      Run expected scan
    --ncpu=<c>       Number of CPUs to use [default: 1]
    --channels=<c1,c2> Channels to use [default: of,sf]

"""

import ROOT as R
import json

from set_limits import *

def create_histfactory(signal_file, prefix, m1, m2, channels):
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
    meas.PrintTree()

    R.RooStats.HistFactory.MakeModelAndMeasurementFast(meas)

def run_limit(sig_file, mass1, mass2, chans, ncpu, asymptotic, expected):
    prefix = "limits/"+sig_file[:-5]+"_{0}_{1}_2".format(mass1, mass2)

    try:
        create_histfactory(sig_file, prefix, mass1, mass2, chans)
    except:
        return None

    if asymptotic:
        res = asymptotic_limit(prefix+"_combined_meas_model.root", expected)
    else:
        res = frequentist_limit(prefix+"_combined_meas_model.root", ncpu, expected)

    return (res.GetExpectedUpperLimit(), res.GetExpectedUpperLimit(-1), res.GetExpectedUpperLimit(+1), res.UpperLimit())

if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)

    if not args['batch']:
        m1 = int(args['<mass1>'])
        m2 = int(args['<mass2>'])
    else:
        with open(args['<mass_file>']) as f:
            masses = json.load(f)
        m1, m2 = masses[int(args['<jobnum>'])]

    sig_file = args['<signal_file>']

    try:
        chans = args['--channels'].split(",")
    except AttributeError:
        chans = ['of', 'sf']

    asym = bool(args['--asymptotic'])

    expected = bool(args['--expected'])

    ncpu = int(args['--ncpu'])

    # reset sys.argv to placate ROOT
    sys.argv = [sys.argv[0], '-b']

    res = run_limit(sig_file, m1, m2, chans, ncpu, asym, expected)

    exp, exp_down, exp_up, obs = res

    if not args['batch']:
        print "95% CL CLs upper limit"
        print "Expected: {0:.2f} ({1:.2f}-{2:.2f})\tObserved: {3:.2f}".format(exp, exp_down, exp_up, obs)
    else:
        with open(args['<output_file>'], 'w') as f:
            if expected:
                f.write("{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4:.2f}\n".format(m1, m2, exp, exp_down, exp_up))
            else:
                f.write("{0}\t{1}\t{2:.2f}\t{3:.2f}\n".format(m1, m2, obs))



