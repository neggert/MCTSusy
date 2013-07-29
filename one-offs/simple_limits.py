import pandas as pd
import selection
from config.parameters import *
import os.path
import ROOT as R
from prep_hists import load_xsec

channels=['of', 'sf']

def get_sms_yield(sms_filename, xsec_filename, hist_filename, mct_cut=120., xsec_multiplier=1.):



    sms_file = pd.HDFStore(sms_filename)
    sms = sms_file['data']
    sel_sms = selection.get_samples(sms)
    mumu_high_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) > 1.)
    mumu_low_eta_sms = sel_sms['opposite_sign_mumu'] & (abs(sms.eta2) < 1.)

    weight = (sel_sms['opposite_sign_ee'].astype(float)*ee_trigger_eff+mumu_high_eta_sms.astype(float)*mumu_high_eta_trigger_eff
                  +mumu_low_eta_sms.astype(float)*mumu_low_eta_trigger_eff + sel_sms['opposite_sign_emu'].astype(float)*emu_trigger_eff)
    weight.name="weight"
    sms = sms.join(weight)

    # check to see if the file exists, since ROOT will happily continue along with a non-existent file
    if not os.path.exists(hist_filename):
        raise IOError(hist_filename+" does not exist.")
    nevents_file = R.TFile(hist_filename)
    nevents_hist = nevents_file.Get("ScanValidator/NEvents_histo")

    xsec_dict = load_xsec(xsec_filename)


    groups = sms.groupby(['mass1', 'mass2'])
    output = []
    for name, data in groups:
        m1, m2 = name
        sel = selection.get_samples(data)

        events_per_point = nevents_hist.GetBinContent(nevents_hist.FindBin(m1, m2))
        try:
            xsec, xsec_err = xsec_dict[float(m1)]
        except KeyError:
            continue

        xsec *= xsec_multiplier

        for ch in channels:
            mass_point = data[sel['sig_'+ch] & (data.mctperp > mct_cut)]

            mass_point.weight *= xsec*lumi/events_per_point
            this_point = {"mass1":m1,
                          "mass2":m2,
                          "channel":ch,
                          "yield": mass_point.weight.sum()}
            output.append(this_point)

    return pd.DataFrame(output)