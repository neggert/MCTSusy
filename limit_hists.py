import sys
import numpy as np
import ROOT as R

from prep_hists import load_xsec

filename = sys.argv[1]

xsec_name = sys.argv[2]

xsecs = load_xsec( xsec_name)

out_filename = sys.argv[3]

data = np.loadtxt(filename)

mass1, mass2 = zip(*data[:,:2])

bin_spacing=25.

m1_low = min(mass1)-bin_spacing/2
m1_high = max(mass1)+bin_spacing/2
m2_low = min(mass2)-bin_spacing/2
m2_high = max(mass2)+bin_spacing/2


nbins_x = int((m1_high-m1_low)/bin_spacing)
nbins_y = int((m2_high-m2_low)/bin_spacing)


# create ROOT file
rfile = R.TFile(out_filename, "RECREATE")
rfile.cd()

# create histograms
xsec_hist = R.TH2D("processXSection", "processXSection", nbins_x, m1_low, m1_high, nbins_y, m2_low, m2_high)
observed_hist = R.TH2D("hObsLimitPlain", "hObsLimitPlain", nbins_x, m1_low, m1_high, nbins_y, m2_low, m2_high)
observed_hist_smooth = R.TH2D("hObsLimitSmooth", "hObsLimitSmooth", nbins_x, m1_low, m1_high, nbins_y, m2_low, m2_high)
expected_hist = R.TH2D("hExpLimitSmooth", "hExpLimitSmooth", nbins_x, m1_low, m1_high, nbins_y, m2_low, m2_high)
expected_m1_hist = R.TH2D("hExpM1LimitSmooth", "hExpM1LimitSmooth", nbins_x, m1_low, m1_high, nbins_y, m2_low, m2_high)
expected_p1_hist = R.TH2D("hExpP1LimitSmooth", "hExpP1LimitSmooth", nbins_x, m1_low, m1_high, nbins_y, m2_low, m2_high)

for d in data:
    m1, m2, exp, exp_m1, exp_p1, obs = d

    xsec_hist.SetBinContent(xsec_hist.FindBin(m1, m2), xsecs[float(m1)][0])
    observed_hist.SetBinContent(observed_hist.FindBin(m1,m2), obs)
    observed_hist_smooth.SetBinContent(observed_hist_smooth.FindBin(m1,m2), obs)
    expected_hist.SetBinContent(expected_hist.FindBin(m1,m2), exp)
    expected_m1_hist.SetBinContent(expected_m1_hist.FindBin(m1,m2), exp_m1)
    expected_p1_hist.SetBinContent(expected_p1_hist.FindBin(m1,m2), exp_p1)

observed_hist_smooth.Smooth(1,"k3a")
expected_hist.Smooth(1,"k3a")
expected_p1_hist.Smooth(1,"k3a")
expected_m1_hist.Smooth(1,"k3a")


xsec_hist.Write()
observed_hist.Write()
observed_hist_smooth.Write()
expected_hist.Write()
expected_m1_hist.Write()
expected_p1_hist.Write()

rfile.Close()



