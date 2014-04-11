#!/usr/bin/env python

"""Prepare histograms

Usage:
    prep_hists2.py <signal_file> <signal_output> <xsec_file> <nevents_file> [-b bins] [-l LOW] [-u high] [-x xsec]

Options:
    -h --help               Show this help information
    -b --bins=<b>              Number of bins [default: 19]
    -l --low=<l>               Lower end of histograms [default: 10]
    -u --high=<h>              Upper end of histograms [default: 200]
    -x --xsec_multiplier=<x>   Cross-section multiplier [default: 1]
"""

from prep_hists import *

from config.data import *


backgrounds = ['of', 'vv', 'wjets', 'z']

def create_template_file(filename="templates.root", bins=19, histrange=(10, 200)):
    """
    Create a ROOT file containing all of the background templates
    """

    wz_sf = (mc.mc_cat=="WZ") & abs(mc.parentParentPdg1).isin([22,23]) & abs(mc.parentParentPdg1).isin([22,23])

    mczz = mc[(mc.mc_cat=='ZZ')]
    mcwz = mc[wz_sf]
    mcvv = mczz.append(mcwz)
    mcz = mc[mc.mc_cat=='DY']
    selvv = selection.get_samples( mcvv, 100.)
    selz = selection.get_samples( mcz, 100.)

    rfile = R.TFile(filename, "RECREATE")
    rfile.cd()

    templates = {}

    of = data[sd['sig_of']]
    templates['of'] = rootutils.create_TH1(of.mctperp, of.weight, "of_template", bins, histrange, True)
    wz_fs = (mc.mc_cat=="WZ") & ~(abs(mc.parentParentPdg1).isin([22,23]) & abs(mc.parentParentPdg1).isin([22,23]))
    fs_mc_yield = mc[smc['sig_mct_low_sf'] & (((mc.mc_cat=="top") & (mc['gen_neutrinos'] >= 2)) | (mc.mc_cat=="WW") | wz_fs | (mc.mc_cat=="HWW") | (mc.mc_cat=="VVV"))].weight.sum()
    templates['of'].Scale(fs_mc_yield)

    wjets = data[sd['wjets_ctrl_sf']]
    templates['wjets'] = rootutils.create_TH1(wjets.mctperp, wjets.weight, "wjets_template", bins, histrange, True)
    
    # systematic on w+jets template
    rhist = R.TH1D("wjets_syst", "wjets_syst", bins, histrange[0], histrange[1])
    for i in xrange(bins):
        # if templates['wjets'].GetBinContent(i+1) > 0: #only do non-zero bins
        rhist.SetBinContent(i+1, 0.3) # 50% systematic
    templates['wjets_syst'] = rhist
    templates['wjets'].Scale(mc[smc['sig_mct_low_sf'] & ((mc['mc_cat']=='fake') | ((mc['mc_cat']=='top') & (mc['gen_neutrinos'] < 2)))].weight.sum())

    vv = mcvv[selvv['sig_sf']]
    templates['vv'] = rootutils.create_TH1(vv.mctperp, vv.weight, "vv_template", bins, histrange, False)


    weights = vv.weight*(1+(vv.mc_cat=="WZ")*0.10)
    templates['vv_syst_WZ_Up'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_WZ_Up", bins, histrange, False)
    weights = vv.weight*(1-(vv.mc_cat=="WZ")*0.10)
    templates['vv_syst_WZ_Down'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_WZ_Down", bins, histrange, False)

    weights = vv.weight*(1+(vv.mc_cat=="ZZ")*0.10)
    templates['vv_syst_ZZ_Up'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_ZZ_Up", bins, histrange, False)
    weights = vv.weight*(1-(vv.mc_cat=="ZZ")*0.10)
    templates['vv_syst_ZZ_Down'] = rootutils.create_TH1(vv.mctperp, weights, "vv_syst_ZZ_Down", bins, histrange, False)

    z = mcz[selz['sig_sf']]
    templates['z'] = rootutils.create_TH1(z.mctperp, z.weight, "z_template", bins, histrange, False)

    # systematic on Z monte carlo
    mc_onz = mc[smc['z_ctrl_sf']]
    data_onz = data[sd['z_ctrl_sf']]

    mc_hist, mc_edges = np.histogram(mc_onz.mctperp, weights=mc_onz.weight, bins=bins, range=histrange, normed=True)
    d_hist, d_edges = np.histogram(data_onz.mctperp, weights=data_onz.weight, bins=bins, range=histrange, normed=True)

    err = np.zeros(bins)
    err[:] = np.max(abs(mc_hist[:10]-d_hist[:10])/d_hist[:10])
    # err[10:] = 0.5

    # make a TH1 out of it
    rhist = R.TH1D("z_syst", "z_syst", bins, histrange[0], histrange[1])
    for i, val in enumerate(err):
        rhist.SetBinContent(i+1, val)

    templates['z_syst'] = rhist


    for k in templates.keys():
        templates[k].Write()

    rfile.Close()

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)

    print(args)

    bins = int(args['--bins'])
    histrange = (float(args['--low']), float(args['--high']))

    create_data_file("data2.root", bins, histrange)
    create_template_file("templates2.root", bins, histrange)
    create_signal_file(args['<signal_file>'], args['<signal_output>'], args['<nevents_file>'], args['<xsec_file>'],
                       float(args['--xsec_multiplier']), bins, histrange)
