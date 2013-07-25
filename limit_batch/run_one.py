#!/usr/bin/env python
"""Run one point

Usage:
    run_one.py batch <batch_file> <batch_num> <label> <output_folder>

Options:

"""

import ROOT as R

def run_one_point(filename, poi_val):
    rfile = R.TFile(filename)
    print filename
    print poi_val

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    sbmodel = ws.obj("ModelConfig")
    sbmodel.GetParametersOfInterest().first().setMax(1000)

    constr = sbmodel.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

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

    toymc = calc.GetTestStatSampler()

    prof_l = R.RooStats.ProfileLikelihoodTestStat(sbmodel.GetPdf())
    prof_l.SetOneSided(True)
    prof_l.SetReuseNLL(True)
    prof_l.SetStrategy(0)

    toymc.SetTestStatistic(prof_l)
    toymc.SetNToys(100) # needed because of https://savannah.cern.ch/bugs/?93360
    toymc.SetGenerateBinned(True)
    toymc.SetUseMultiGen(True)

    hypo.RunOnePoint(poi_val)

    res = hypo.GetInterval()

    return res

if __name__ == '__main__':
    from docopt import docopt
    import re
    args = docopt(__doc__)

    print args

    if args['batch']:
        with open(args['<batch_file>']) as bf:
            model_file, poi = bf.readlines()[int(args['<batch_num>'])/50].split()
        m1, m2 = map(int, re.search("_(\d+)_(\d+)_", model_file).groups())
        output_file = args['<output_folder>']+"/{0}_{1}_{2}_exp{3}.root".format(args["<label>"], m1, m2, int(args['<batch_num>']))
    else :
        model_file = args['<model_file>']
        poi = args['<poi_val>']
        output_file = args['<output_file>']

    res = run_one_point(model_file, float(poi))

    f = R.TFile(output_file, "RECREATE")
    f.cd()
    res.Write()
    f.Close()

