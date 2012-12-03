import selection
import rootutils

from collections import defaultdict

# <codecell>

import ROOT as r
r.gSystem.SetDynamicPath(r.gSystem.Getenv("LD_LIBRARY_PATH"))
r.gSystem.Load("libRooFitCore")
r.gSystem.Load("libRooFit")
r.gSystem.Load("libRooStats")


class StatModel(object):
    """Statistical model for MCT analysis"""
    def __init__(self, data, mc, signal, mct_cut):
        """
        Initialize model

        Arguments:
        data - dataframe containing observations
        mc - dataframe containing simulation
        signal - dataframe containing the signal
        mct_cut - float defining MCT cut for signal region
        """

        super(StatModel, self).__init__()

        self.mct_cut = mct_cut
        self.channels = ['_of', '_sf']
        self.roofit_channel = r.RooCategory("channel", "channel")
        for c in self.channels:
            self.roofit_channel.defineType(c)
        self.backgrounds = ['top', 'wjets', 'vv', 'z']

        # set up the data samples
        self.data = data
        self.mcvv = mc[(mc.mc_cat=='WV') | (mc.mc_cat=='ZZ')]
        self.mcz = mc[mc.mc_cat=='DY']

        # set up the selections
        self.sel = selection.get_samples(self.data, mct_cut)
        self.selvv = selection.get_samples(self.mcvv, mct_cut)
        self.selz = selection.get_samples(self.mcz, mct_cut)

        # set up RooFit variables
        self.mct = r.RooRealVar("mct", "mct", 5., 100.)  # TODO: fix upper limit
        self.w = r.RooRealVar("w", "w", 0., 10.)
        self.evt_yield = r.RooRealVar("yield", "yield", 0., 1000.)

        # number of events in the sideband regions
        # one copy for each channel
        self.n_sb = defaultdict(dict)
        self.n_ctrl_sb = defaultdict(dict)
        self.n_ctrl_sig = defaultdict(dict)
        for bkg in self.backgrounds:
            for channel in self.channels:
                # n_sb comes from the normalization fit in the sideband
                self.n_sb[bkg][channel] =       r.RooRealVar("n_sb_{0}_{1}".format(bkg, channel),
                                                             "n_sb_{0}_{1}".format(bkg, channel),
                                                             1000, -10000., 10000.)
                # This comes from the sideband in the control sample
                # It's treated as having no uncertainty
                self.n_ctrl_sb[bkg][channel] =  r.RooRealVar("n_ctrl_sb_{0}_{1}".format(bkg, channel),
                                                             "n_ctrl_sb_{0}_{1}".format(bkg, channel),
                                                             1000, -10000., 10000.)
                # This comes from the high-mct signal region in the control sample
                # It's treated as a sample from a Poisson-distributed quantity with unknown mean
                self.n_ctrl_sig[bkg][channel] = r.RooRealVar("n_ctrl_sig_{0}_{1}".format(bkg, channel),
                                                             "n_ctrl_sig_{0}_{1}".format(bkg, channel),
                                                             1000, -10000., 10000.)

        # set up dataset
        ds_ch = {}

        for channel in self.channels:
            ds_ch[channel] = rootutils.create_roodataset(self.data[self.sel['sig_mct_low'+channel]].mctperp,
                                                self.data[self.sel['sig_mct_low'+channel]].weight,
                                                self.mct, self.w, title="data")

        self.dataset = r.RooDataSet("dataset", "dataset", r.RooArgSet(self.mct, self.w), r.RooFit.Index(self.roofit_channel),
                                    r.RooFit.Import("_of", ds_ch['_of']),
                                    r.RooFit.Import("_sf", ds_ch['_sf'])
                                   )

    def build_sideband_pdf(self):
        """
        Build the likelihood for the normalization fit in the sideband region

        Arguments:
        mct_range - range defining the MCT sideband

        returns: RooAbsPdf - likelihood function for sideband normalizations
        """

        # Some of these we won't need again, but we make them into class members
        # to keep them from getting garbage-collected, which causes segfaults.
        # At least, I think that's what's going on.
        self.models_by_channel = {}
        self.shapehistTH1 = defaultdict(dict)
        self.shapehist = defaultdict(dict)
        self.shape = defaultdict(dict)
        self.shape_ext = defaultdict(dict)

        # do this for both channels:
        sim_pdf = r.RooSimultaneous("sideband_normalization", "sideband_normalization", self.roofit_channel)

        for channel in self.channels:
            sum_pdfs = r.RooArgList()
            # convert data to TH1s
            self.shapehistTH1['top'][channel] = rootutils.create_TH1(self.data[self.sel['top_mct_low'+channel]].mctperp, self.data[self.sel['top_mct_low'+channel]].weight, "top")
            self.shapehistTH1['wjets'][channel] = rootutils.create_TH1(self.data[self.sel['wjets_mct_low'+channel]].mctperp, self.data[self.sel['wjets_mct_low'+channel]].weight, "wjets")
            self.shapehistTH1['vv'][channel] = rootutils.create_TH1(self.mcvv[self.selvv['sig_mct_low'+channel]].mctperp, self.mcvv[self.selvv['sig_mct_low'+channel]].weight, "vv")
            self.shapehistTH1['z'][channel] = rootutils.create_TH1(self.mcz[self.selz['sig_mct_low'+channel]].mctperp, self.mcz[self.selz['sig_mct_low'+channel]].weight, "Z")

            for bkg in self.backgrounds:
                # convert TH1 to a RooDataHist
                self.shapehist[bkg][channel] = r.RooDataHist("hist_{}_{}".format(bkg, channel),
                                                "hist_{}_{}".format(bkg, channel),
                                                r.RooArgList(self.mct), self.shapehistTH1[bkg][channel])
                # make a PDF out of the RooDataHist
                self.shape[bkg][channel] = r.RooHistPdf("shape_{}_{}".format(bkg, channel),
                                                        "shape_{}_{}".format(bkg, channel),
                                                         r.RooArgSet(self.mct), self.shapehist[bkg][channel])
                # turn the PDF into an extended PDF
                self.shape_ext[bkg][channel] = r.RooExtendPdf("shape_ext_{}_{}".format(bkg, channel),
                                                              "shape_ext_{}_{}".format(bkg, channel),
                                                              self.shape[bkg][channel],
                                                              self.n_sb[bkg][channel])
                sum_pdfs.add(self.shape_ext[bkg][channel])


            self.models_by_channel[channel] = r.RooAddPdf("model"+channel, "model"+channel, sum_pdfs)
            sim_pdf.addPdf(self.models_by_channel[channel], channel)

        return sim_pdf

    def fit_sb_normalization(self):
        """This is mostly for testing purposes"""
        model = self.build_sideband_likelihood()

        return model.fitTo(self.dataset, r.RooFit.Save())




    def build_signal_region_prediction(self):
        """
        Build the predictions for the number of events in the signal region
        """
        self.ctrl_ratio = defaultdict(dict)
        self.n_sig = defaultdict(dict)

        for channel in self.channels:
            for bkg in self.backgrounds:
                self.ctrl_ratio[bkg][channel] = r.RooFormulaVar("ctrl_ratio_{}_{}".format(bkg, channel),
                                                           "ctrl_ratio_{}_{}".format(bkg, channel),
                                                           "n_sb_{0}_{1}/n_ctrl_sb_{0}_{1}".format(bkg, channel),
                                                           r.RooArgList(self.n_sb[bkg][channel], self.n_ctrl_sb[bkg][channel]))
                self.n_sig[bkg][channel] = r.RooFormulaVar("n_sig_{}_{}".format(bkg, channel),
                                                      "n_sig_{}_{}".format(bkg, channel),
                                                      "ctrl_ratio_{0}_{1}*n_ctrl_sig_{0}_{1}".format(bkg, channel),
                                                      r.RooArgList(self.ctrl_ratio[bkg][channel], self.n_ctrl_sig[bkg][channel])
                                                      )

    def build_signal_region_pdf(self):
        """
        Build the pdf for the signal region
        """
        self.build_signal_region_prediction()

        sig_pdf = r.RooSimultaneous("sig_pdf", "sig_pdf", self.roofit_channel)
        self.sig_pdfs_by_channel = {}
        self.sig_bkg = {}
        for channel in self.channels:
            bkg_formula = ""
            bkg_contribs = r.RooArgList()
            for bkg in self.backgrounds:
                bkg_formula += "n_sig_{}_{}+".format(bkg, channel)
                bkg_contribs.add(self.n_sig[bkg][channel])
            bkg_formula = bkg_formula[:-1]  # chop off the trailing +
            self.sig_bkg[channel] = r.RooFormulaVar("sig_bkg"+channel, "sig_bkg"+channel, bkg_formula, bkg_contribs)
            self.sig_pdfs_by_channel[channel] = r.RooPoisson("signal_pmf"+channel, "signal_pmf"+channel, self.evt_yield, self.sig_bkg[channel])
            sig_pdf.addPdf(self.sig_pdfs_by_channel[channel], channel)

        return sig_pdf

