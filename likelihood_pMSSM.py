#!/usr/bin/env python

import ROOT as R
import json
import re

R.SetMemoryPolicy(R.kMemoryStrict)

with open("limits/nevents_pMSSM.json") as f:
    nevents = json.load(f)


def get_likelihood(model_file, mu, mu_fixed=True):
    rfile = R.TFile(model_file)

    ws = rfile.Get("combined")

    data = ws.data("obsData")

    model = ws.obj("ModelConfig")

    constr = model.GetNuisanceParameters()
    R.RooStats.RemoveConstantParameters(constr)

    model.GetParametersOfInterest().first().setVal(mu)
    if mu_fixed:
        model.GetParametersOfInterest().first().setConstant()


    # run the fit
    res = model.GetPdf().fitTo(data, R.RooFit.Constrain(constr), R.RooFit.Save(), R.RooFit.PrintLevel(0))
   
    ll = -res.minNll()

    res.IsA().Destructor(res)
    constr.IsA().Destructor(constr)

    rfile.Close()

    model.IsA().Destructor(model)
    data.IsA().Destructor(data)

    #ws.IsA().Destructor(ws)

    #rfile.IsA().Destructor(rfile)

    return ll

def process_model(model_file):
    point_num = re.search("_(\d*?)_", model_file).groups()[0]
    n = nevents[point_num]
    maxllhd = get_likelihood(model_file, 1., False)
    maxllhd_000 = get_likelihood(model_file, 0.)
    maxllhd_050 = get_likelihood(model_file, 0.5)
    maxllhd_100 = get_likelihood(model_file, 1.0)
    maxllhd_150 = get_likelihood(model_file, 1.5)

    return point_num, n, maxllhd_000, maxllhd_050, maxllhd_100, maxllhd_150, maxllhd

if __name__ == '__main__':
    import sys
    files = sys.argv[2:]
    with open(sys.argv[1], 'w') as f:
        f.write("ID\tN_TOT\tmaxllhd_000\tmaxllhd_050\tmaxllhd_100\tmaxllhd_150\tmaxllhd\n")
        for model_file in files:
            data = process_model(model_file)
            f.write("\t".join(map(str, data)) + "\n")
            R.gROOT.Clear()
            R.gROOT.DeleteAll()


