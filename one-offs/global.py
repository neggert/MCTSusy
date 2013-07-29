import ROOT
import numpy as np

def print_arg_list(l):
    for i in xrange(l.getSize()):
        l.at(i).Print()

f = ROOT.TFile("limits/sig_chi_600_200_combined_meas_model.root")

ws = f.Get("combined")

model = ws.obj("ModelConfig")
data = ws.data("obsData")

pdf = model.GetPdf()
go_1 = ws.obj("nom_gamma_ww_syst_sf_bin_11")
go_1.setVal(0.5)

pdf.fitTo(data, ROOT.RooFit.PrintLevel(0))

# dummy = ROOT.RooStats.MinNLLTestStat(model.GetPdf())
# mc = ROOT.RooStats.ToyMCSampler(dummy, 1)
# mc.SetPdf(model.GetPdf())
# mc.SetObservables(model.GetObservables())
# mc.SetGlobalObservables(model.GetGlobalObservables())
# mc.SetParametersForTestStat(model.GetParametersOfInterest())
# 
# go_2 = ws.obj("nom_gamma_ww_syst_sf_bin_10")
# 
# results = []
# 
# for i in xrange(1000):
#     mc.GenerateGlobalObservables(pdf)
#     results.append((go_1.getVal(), go_2.getVal()))
# 
# a = np.asarray(results)
# 
# print np.cov(a, rowvar=0)

