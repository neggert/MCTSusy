import sys
import numpy as np
import ROOT as R

def load_xsec(filename):
    """
    load cross-sections for a file. Returns a dictionary describing cross-sections on their uncertainties.
    The dictionary is indexed by the particle mass, and each element is a tuple containing (xsec, uncertainty)
    """
    f = open(filename)
    xsec_dict = {}
    for line in f:
        mass, xsec, err = line.split()
        xsec_dict[float(mass)] = (float(xsec), float(err))
    f.close()
    return xsec_dict

filename = sys.argv[1]

xsec_name = sys.argv[2]

xsecs = load_xsec( xsec_name)

out_filename = sys.argv[3]

data = np.loadtxt(filename, delimiter=",")

mass1, mass2 = zip(*data[:,:2])

bin_spacing=25.

m1_low = 100-bin_spacing/2
m1_high = 775+bin_spacing/2
m2_low = 0-bin_spacing/2
m2_high = 625+bin_spacing/2


nbins_x = int((m1_high-m1_low)/bin_spacing)
nbins_y = int((m2_high-m2_low)/bin_spacing)


# create ROOT file
# rfile = R.TFile(out_filename, "RECREATE")
# rfile.cd()

# create histograms
expected_hist = R.TH2D("hExpLimitSmooth", "hExpLimitSmooth", nbins_x, m1_low, m1_high, nbins_y, m2_low, m2_high)

the_data = {}
for d in data:
    m1, m2, exp = d
    the_data[(m1,m2)] = exp

import json
with open("chi_masses.json") as f:
    masses = json.load(f)

for m1, m2 in masses:

    if m1-m2 <= 50:
        continue

    try:
        exp = the_data[(m1, m2)]
    except KeyError:
        continue

    expected_hist.SetBinContent(expected_hist.FindBin(m1,m2), exp)


expected_hist.Smooth(1,"k3a")

old_file = R.TFile("chipm_limits_20130503.root")
fullHist = old_file.Get("hExpLimitSmooth")

fullHist.Divide(expected_hist)

fullHist.GetXaxis().SetTitle("mChargino")
fullHist.GetYaxis().SetTitle("mLSP")
fullHist.GetZaxis().SetRangeUser(0., 1.)

c1 = R.TCanvas()
fullHist.Draw("colz")
import IPython
IPython.embed()
c1.SaveAs("chipm_limit_ratio.png")

# expected_hist.Write()

# rfile.Close()



