#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    set_limits.py <signal_file> <mass1> <mass2> [-hac] [--ncpu=<c>] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    -a --asymptotic  Run asymptotic limits
    -c --coarse      Run coarse scan
    --ncpu=<c>       Number of CPUs to use [default: 1]
    --channels=<c1,c2> Channels to use [default: of,sf]

"""

backgrounds = ['top', 'vv', 'wjets', 'z']

import ROOT as R
from collections import defaultdict

def create_histfactory(signal_file, prefix, m1, m2, channels):
    meas = R.RooStats.HistFactory.Measurement("meas", "meas")

    meas.SetOutputFilePrefix(prefix)
    meas.SetPOI("sig_strength")
    meas.AddConstantParam("n_of_top")
    meas.AddConstantParam("n_sf_top")

    temp_file = R.TFile("templates.root")

    meas.SetLumi(1.0)
    meas.SetLumiRelErr(0.04)
    meas.SetExportOnly(True)

    channel_confs = {}
    samples = defaultdict(dict)

    for ch in channels:
        channel_confs[ch] = R.RooStats.HistFactory.Channel(ch)
        channel_confs[ch].SetData("data_"+ch, "data.root")
        channel_confs[ch].SetStatErrorConfig(0.05, "Poisson")

        # signal sample
        signal = R.RooStats.HistFactory.Sample("signal_"+ch, "sms_template_{}_{}_{}".format(ch, m1, m2), signal_file)
        signal.SetNormalizeByTheory(True)
        signal.AddNormFactor("sig_strength", 0.1, 0., 10)
        signal.ActivateStatError()
        signal.AddOverallSys("trigger", 0.95, 1.05)
        signal.AddOverallSys("id_and_selection", 0.98, 1.02)
        signal.AddOverallSys("b_veto", 0.94, 1.06)
        signal.AddHistoSys("jes", "sms_template_jes_down_{}_{}_{}".format(ch, m1, m2), signal_file, "",
                           "sms_template_jes_up_{}_{}_{}".format(ch, m1, m2), signal_file, "")
        channel_confs[ch].AddSample(signal)

        # add the background samples
        for bkg in backgrounds:
            template = R.RooStats.HistFactory.Sample("{}_{}".format(bkg, ch), "{}_template_{}".format(bkg, ch), "templates.root")
            template.SetNormalizeByTheory(False)
            template.ActivateStatError()
            if bkg == 'top':
                ntop_pred = temp_file.Get("ntop_"+ch)[0]
                template.AddNormFactor("n_{}_{}".format(ch, bkg), ntop_pred, 0, 2*ntop_pred, True)
                template.AddOverallSys("top_norm_"+ch, 0.88, 1.12)
            else :
                template.AddNormFactor("n_{}_{}".format(ch, bkg), 2000, 0, 5000)


            samples[ch][bkg] = template
            channel_confs[ch].AddSample(samples[ch][bkg])

        meas.AddChannel(channel_confs[ch])

    meas.CollectHistograms()
    meas.PrintTree()

    R.RooStats.HistFactory.MakeModelAndMeasurementFast(meas)

def asymptotic_limit(filename, coarse):
    rfile = R.TFile(filename)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    sbmodel = ws.obj("ModelConfig")

    constr = sbmodel.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    # run the initial fit
    sbmodel.GetPdf().fitTo(data, R.RooFit.Constrain(constr))

    sbmodel.SetSnapshot(sbmodel.GetParametersOfInterest())
    sbmodel.Print()

    bmodel = sbmodel.Clone()

    bmodel.GetParametersOfInterest().first().setVal(0.)
    bmodel.SetName("ModelConfig_bonly")
    bmodel.SetSnapshot(bmodel.GetParametersOfInterest())

    R.RooStats.AsymptoticCalculator.SetPrintLevel(0)
    calc = R.RooStats.AsymptoticCalculator(data, bmodel, sbmodel)
    calc.SetOneSided(True)

    hypo = R.RooStats.HypoTestInverter(calc)
    hypo.SetConfidenceLevel(0.95)
    hypo.UseCLs(True)
    hypo.SetVerbose(True)

    if coarse:
        hypo.SetFixedScan(10, 0, 10)

    res = hypo.GetInterval()

    return res

def frequentist_limit(filename, ncpu, coarse):
    rfile = R.TFile(filename)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    sbmodel = ws.obj("ModelConfig")

    constr = sbmodel.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    # run the initial fit
    sbmodel.GetPdf().fitTo(data, R.RooFit.Constrain(constr))

    sbmodel.SetSnapshot(sbmodel.GetParametersOfInterest())

    bmodel = sbmodel.Clone()

    bmodel.GetParametersOfInterest().first().setVal(0.)
    bmodel.SetName("ModelConfig_bonly")
    bmodel.SetSnapshot(bmodel.GetParametersOfInterest())

    calc = R.RooStats.FrequentistCalculator(data, bmodel, sbmodel)

    bmodel.Print()
    sbmodel.Print()

    hypo = R.RooStats.HypoTestInverter(calc)
    hypo.UseCLs(True)
    hypo.SetVerbose(True)

    if coarse:
        hypo.SetFixedScan(10, 0, 10)

    toymc = calc.GetTestStatSampler()

    prof_l = R.RooStats.ProfileLikelihoodTestStat(sbmodel.GetPdf())
    prof_l.SetOneSided(True)

    toymc.SetTestStatistic(prof_l)

    if ncpu > 1:
        pc = R.RooStats.ProofConfig(ws, ncpu, "")
        toymc.SetProofConfig(pc)

    res = hypo.GetInterval()

    return res

if __name__ == '__main__':
    from docopt import docopt

    args = docopt(__doc__)

    print args

    m1 = int(args['<mass1>'])
    m2 = int(args['<mass2>'])

    sig_file = args['<signal_file>']

    prefix = "limits/"+sig_file[:-5]+"_{}_{}".format(m1, m2)
    try:
        chans = args['--channels'].split(",")
    except AttributeError:
        chans = ['of', 'sf']


    create_histfactory(sig_file, prefix, m1, m2, chans)

    asym = bool(args['--asymptotic'])

    coarse = bool(args['--coarse'])

    if asym:
        res = asymptotic_limit(prefix+"_combined_meas_model.root", coarse)
    else:
        res = frequentist_limit(prefix+"_combined_meas_model.root", int(args['--ncpu']), coarse)

    exp = res.GetExpectedUpperLimit()
    obs = res.UpperLimit()

    print "95% CL CLs upper limit"
    print "Expected: {:.2f}\tObserved: {:.2f}".format(exp, obs)


