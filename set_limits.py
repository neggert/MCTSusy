#! /usr/bin/env python

"""Set limits on the specified model

Usage:
    set_limits.py <signal_file> <mass1> <mass2> [-hac] [--ncpu=<c>] [--channels=<c1,c2>]
    set_limits.py batch <signal_file> <mass_file> <jobnum> <output_file> [-ah] [--channels=<c1,c2>]

Options:
    -h --help        Show this screen.
    -a --asymptotic  Run asymptotic limits
    -c --coarse      Run coarse scan
    --ncpu=<c>       Number of CPUs to use [default: 1]
    --channels=<c1,c2> Channels to use [default: of,sf]

"""
backgrounds = ['top', 'vv', 'wjets', 'z']

import ROOT as R
import json
from collections import defaultdict

def create_histfactory(signal_file, prefix, m1, m2, channels, data_file_name="data.root", data_prefix="data"):
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
            if bkg == 'banana':
                ntop_pred = temp_file.Get("ntop_"+ch)[0]
                template.AddNormFactor("n_{0}_{1}".format(ch, bkg), ntop_pred, 0, 2*ntop_pred, True)
                template.AddOverallSys("top_norm_"+ch, 0.88, 1.12)
            else :
                template.AddNormFactor("n_{0}_{1}".format(ch, bkg), 2000, 0, 10000)

            if bkg == 'z':
                template.AddShapeSys("z_syst_"+ch, 0, "z_syst", "templates.root")
            if bkg == 'wjets':
                template.AddShapeSys("wjets_syst_"+ch, 0, "wjets_syst_"+ch, "templates.root")


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

    poi_hat = sbmodel.GetParametersOfInterest().first().getVal()
    poi_hat_err = sbmodel.GetParametersOfInterest().first().getError()

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
    hypo.SetFixedScan(50, poi_hat, poi_hat+4*poi_hat_err, True)

    res = hypo.GetInterval()

    return res

def frequentist_limit(filename, ncpu, coarse):
    # First run the asymptotic limit to get a rough idea
    rfile = R.TFile(filename)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    sbmodel = ws.obj("ModelConfig")

    constr = sbmodel.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    # run the initial fit
    sbmodel.GetPdf().fitTo(data, R.RooFit.Constrain(constr))

    poi_hat = sbmodel.GetParametersOfInterest().first().getVal()
    poi_hat_err = sbmodel.GetParametersOfInterest().first().getError()

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
    else:
        hypo.SetFixedScan(50, poi_hat, poi_hat+4*poi_hat_err, True)

    toymc = calc.GetTestStatSampler()

    prof_l = R.RooStats.ProfileLikelihoodTestStat(sbmodel.GetPdf())
    prof_l.SetOneSided(True)

    toymc.SetTestStatistic(prof_l)
    toymc.SetMaxToys(500) # needed because of https://savannah.cern.ch/bugs/?93360

    if ncpu > 1:
        pc = R.RooStats.ProofConfig(ws, ncpu, "")
        toymc.SetProofConfig(pc)

    res = hypo.GetInterval()

    return res

def run_limit(sig_file, mass1, mass2, chans, ncpu, asymptotic, coarse):
    prefix = "limits/"+sig_file[:-5]+"_{0}_{1}".format(mass1, mass2)

    # try:
    create_histfactory(sig_file, prefix, mass1, mass2, chans)
    # except:
        # return None

    if asymptotic:
        res = asymptotic_limit(prefix+"_combined_meas_model.root", coarse)
    else:
        res = frequentist_limit(prefix+"_combined_meas_model.root", ncpu, coarse)

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
        m1, m2 = [int(m) for m in masses[int(args['<jobnum>'])]]

    sig_file = args['<signal_file>']

    try:
        chans = args['--channels'].split(",")
    except AttributeError:
        chans = ['of', 'sf']

    asym = bool(args['--asymptotic'])

    coarse = bool(args['--coarse'])

    ncpu = int(args['--ncpu'])

    res = run_limit(sig_file, m1, m2, chans, ncpu, asym, coarse)

    exp, exp_down, exp_up, obs = res

    if not args['batch']:
        print "95% CL CLs upper limit"
        print "Expected: {0:.2f} ({1:.2f}-{2:.2f})\tObserved: {3:.2f}".format(exp, exp_down, exp_up, obs)
    else:
        with open(args['<output_file>'], 'w') as f:
            f.write("{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4:.2f}\t{5:.2f}\n".format(m1, m2, exp, exp_down, exp_up, obs))


