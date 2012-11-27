"""
functions to calculate and write SMS datacards. The major that will be used is
write_sms_dc
"""

from math import sqrt
import sys
sys.path.append("../")
from selection import *
from parameters import *
from data import *
from os.path import exists

import ROOT

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


def write_sms_dc(sms, xsec_filename, hist_filename, outfilename, channels, xsec_scale=1., fixed_xsec=None):
    """
    Calculate and write the SMS datacards

    Keyword Arguments:
    sms -- pandas dataframe containing all of the sms data with no cuts applied
    xsec_filename -- File containing the cross sections as a function of mass. Passed to load_xsec
    histfilename -- ROOT file containing a TH2D (hist/NEvents_histo) giving the number of events generated at each mass point
    outfilename -- Base filename for output cards. The script will produce output files named outfilename_$m1_$m2.txt
    channels -- List or tuple of channels to use. Usually come combination of "_of" and "_sf"
    xsec_scale -- Scale factor for cross sections. Values read in from the file are mulitplied by this. (default 1.0)
    fixed_xsec -- If set, use this number as the cross-section at all mass points instead of reading them from a file
    """

    if fixed_xsec is None:
        xsec_dict = load_xsec(xsec_filename)

    # check to see if the file exists, since ROOT will happily continue along with a non-existent file
    if not exists(hist_filename):
        raise IOError(hist_filename+" does not exist.")
    rfile = ROOT.TFile(hist_filename)
    nevents_hist = rfile.Get("hist/NEvents_histo")

    groups = sms.groupby(['mass1', 'mass2'])

    for name, data in groups:
        m1, m2 = name

        events_per_point = nevents_hist.GetBinContent(nevents_hist.FindBin(m1, m2))

        if fixed_xsec is None:
            try:
                xsec, xsec_err = xsec_dict[float(m1)]
            except KeyError:
                continue
        else :
            xsec, xsec_err = fixed_xsec, 0.
        xsec *= xsec_scale

        f = open(outfilename+"_"+str(int(m1))+"_"+str(int(m2))+".txt", 'w')

        line1 = [str(i) for i in ("xsec", m1, m2, xsec)]
        line1 = " ".join(line1)
        f.write(line1)
        f.write("\n")

        sel = get_samples(data)
        for i, ch in enumerate(channels):
            num_ee = data[sel['sig_mct_high'+ch] & sel['sig_ee']].event.count()
            num_mumu_high_eta = data[sel['sig_mct_high'+ch] & sel['sig_mumu'] & (abs(data.eta2) > 1.)].event.count()
            num_mumu_low_eta = data[sel['sig_mct_high'+ch] & sel['sig_mumu'] & (abs(data.eta2) < 1.)].event.count()
            num_emu = data[sel['sig_mct_high'+ch] & sel['sig_emu']].event.count()

            # scale by trigger efficiency
            num_ee_cor = num_ee*ee_trigger_eff
            num_mumu_high_cor = num_mumu_high_eta*mumu_high_eta_trigger_eff
            num_mumu_low_cor = num_mumu_low_eta*mumu_low_eta_trigger_eff
            num_emu_cor = num_emu*emu_trigger_eff
            num_mumu_cor = num_mumu_high_cor+num_mumu_low_cor

            # eventually need actual number of events per point
            eff = (num_ee_cor+num_mumu_cor+num_emu_cor)/events_per_point

            my_yield = xsec*eff

            # do it scaled up and down
            num_ee_scaleup = data[sel['sig_mct_high_scaleup'+ch] & sel['ee']].event.count()
            num_mumu_high_eta_scaleup = data[sel['sig_mct_high_scaleup'+ch] & sel['mumu'] & (abs(data.eta2) > 1.)].event.count()
            num_mumu_low_eta_scaleup = data[sel['sig_mct_high_scaleup'+ch] & sel['mumu'] & (abs(data.eta2) < 1.)].event.count()
            num_emu_scaleup = data[sel['sig_mct_high_scaleup'+ch] & sel['emu']].event.count()

            num_ee_scaledown = data[sel['sig_mct_high_scaledown'+ch] & sel['ee']].event.count()
            num_mumu_high_eta_scaledown = data[sel['sig_mct_high_scaledown'+ch] & sel['mumu'] & (abs(data.eta2) > 1.)].event.count()
            num_mumu_low_eta_scaledown = data[sel['sig_mct_high_scaledown'+ch] & sel['mumu'] & (abs(data.eta2) < 1.)].event.count()
            num_emu_scaledown = data[sel['sig_mct_high_scaledown'+ch] & sel['emu']].event.count()

            yield_up = get_yield(num_ee_scaleup, num_mumu_high_eta_scaleup, num_mumu_low_eta_scaleup, num_emu_scaleup, events_per_point, xsec)
            yield_down = get_yield(num_ee_scaledown, num_mumu_high_eta_scaledown, num_mumu_low_eta_scaledown, num_emu_scaledown, events_per_point, xsec)

            jes_unc = 1.*abs(yield_up-yield_down)/2

            theory_unc = xsec_err*eff

            trigger_unc = xsec/events_per_point*sqrt( (num_ee_cor*ee_trigger_eff_frac_unc)**2
                                                     +(num_mumu_cor*mumu_trigger_eff_frac_unc)**2
                                                     +(num_emu_cor*emu_trigger_eff_frac_unc)**2)

            ele_selection_unc = xsec*(num_ee_cor+num_emu_cor)/events_per_point*ele_selection_eff_frac_unc
            mu_selection_unc = xsec*(num_mumu_cor+num_emu_cor)/events_per_point*mu_selection_eff_frac_unc
            tau_selection_unc = 0.

            uncorrelated_unc = (sqrt(num_ee)*ee_trigger_eff+sqrt(num_mumu_low_eta)*mumu_low_eta_trigger_eff+
                                sqrt(num_mumu_high_eta)*mumu_high_eta_trigger_eff+sqrt(num_emu)*emu_trigger_eff)/events_per_point*xsec
            assert(uncorrelated_unc <= my_yield)

            total_unc = sqrt(uncorrelated_unc**2
                            +ele_selection_unc**2
                            +mu_selection_unc**2
                            +trigger_unc**2
                            +jes_unc**2
                            +theory_unc**2)


            # scan <m1> <m2> <channel id> <yield in pb> <total yield uncertainty> <uncorrelated uncertainty> <muon eff uncertainty> <electron eff uncertainty> <tau eff uncertainty> <trigger eff uncertainty> <JES uncertainty> <theory uncertainties>
            line2 = [str(i) for i in ("scan", m1, m2, i, my_yield, total_unc, uncorrelated_unc, mu_selection_unc, ele_selection_unc, tau_selection_unc, trigger_unc, jes_unc, theory_unc)]
            line2 = " ".join(line2)
            f.write(line2)
            f.write("\n")
        f.close()

def get_yield(num_ee, num_mumu_low_eta, num_mumu_high_eta, num_emu, events_per_point, xsec):
    """
    Get MC yield in pb

    Arguments:
    num_ee -- number of ee events passing cuts
    num_mumu_low_eta -- number of mumu events at low eta (note: eta affects the trigger efficiency that gets applied)
    num_mumu_high_eta -- number of mumu events at high eta
    num_emu -- number of emu events
    events_per_point -- the number of mc events generated at this mass point. Needed as the denominator when calculating efficiency
    xsec -- cross-section in pb for this mass point
    """
    # scale by trigger efficiency
    num_ee_cor = num_ee*ee_trigger_eff
    num_mumu_high_cor = num_mumu_high_eta*mumu_high_eta_trigger_eff
    num_mumu_low_cor = num_mumu_low_eta*mumu_low_eta_trigger_eff
    num_emu_cor = num_emu*emu_trigger_eff
    num_mumu_cor = num_mumu_high_cor+num_mumu_low_cor

    eff = (num_ee_cor+num_mumu_cor+num_emu_cor)/events_per_point

    return xsec*eff

if __name__ == '__main__':
    # sms = HDFStore("Data/sms_chi.hdf5")
    # data = sms['smsTChipmSlepSnu']

    # write_sms_dc(data, "8TeVc1c1.xsec", "datacards/TChipmSlepSnu/TChipmSlepSnu")

    write_sms_dc(slep, "limits/8TeVeLeL.xsec", "limits/histo_slepslep.root", "limits/datacards/TSlepSlep/TSlepSlep", ["_sf",], xsec_scale=2.)
