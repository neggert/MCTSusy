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

    mcvv = mc[(mc.mc_cat=='ZZ')]
    mcz = mc[mc.mc_cat=='DY']
    selvv = selection.get_samples( mcvv, 100.)
    selz = selection.get_samples( mcz, 100.)

    rfile = R.TFile(filename, "RECREATE")
    rfile.cd()

    templates = {}

    of = data[sd['sig_of']]
    templates['of'] = rootutils.create_TH1(of.mctperp, of.weight, "of_template", bins, histrange, True)

    wjets = data[sd['wjets_ctrl_sf']]
    templates['wjets'] = rootutils.create_TH1(wjets.mctperp, wjets.weight, "wjets_template", bins, histrange, True)

    vv = mcvv[selvv['sig_sf']]
    templates['vv'] = rootutils.create_TH1(vv.mctperp, vv.weight, "vv_template", bins, histrange, True)

    z = mcz[selz['sig_sf']]
    templates['z'] = rootutils.create_TH1(z.mctperp, z.weight, "z_template", bins, histrange, True)

    # systematic on Z monte carlo
    mc_onz = mc[smc['z_ctrl_sf']]
    data_onz = data[sd['z_ctrl_sf']]

    mc_hist, mc_edges = np.histogram(mc_onz.mctperp, weights=mc_onz.weight, bins=bins, range=histrange, normed=True)
    d_hist, d_edges = np.histogram(data_onz.mctperp, weights=data_onz.weight, bins=bins, range=histrange, normed=True)

    err = abs(mc_hist[:10]-d_hist[:10])

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

    create_data_file("data.root", bins, histrange)
    create_template_file("templates2.root", bins, histrange)
    create_signal_file(args['<signal_file>'], args['<signal_output>'], args['<nevents_file>'], args['<xsec_file>'],
                       float(args['--xsec_multiplier']), bins, histrange)
