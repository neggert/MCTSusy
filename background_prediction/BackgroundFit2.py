
# Load libraries and the data

# <codecell>

import sys
sys.path.append("../")
sys.path.append("../import/")

from pandas import *
from math import *
import scipy.stats
from CMSPyLibs.plot import *
from selection import *
from rootutils import *
from collections import defaultdict

import json

# <codecell>

import ROOT as r
r.gSystem.SetDynamicPath(r.gSystem.Getenv("LD_LIBRARY_PATH"));
r.gSystem.Load("libRooFitCore")
r.gSystem.Load("libRooFit")
r.gSystem.Load("libRooStats")

# <codecell>

def do_bkg_fit2(data, mc, mctcut=100., flavor='sf', plot=False) :

      wz_sf = (mc.mc_cat=="WZ") & abs(mc.parentParentPdg1).isin([22,23]) & abs(mc.parentParentPdg1).isin([22,23])

      mczz = mc[(mc.mc_cat=='ZZ') | wz_sf]
      mcz = mc[mc.mc_cat=='DY']

      # <markdowncell>

      # Run selections on the data. This creates all the control and signal samples. Results are stored in sel. More detailed explanations of the different selections are in the docstring for get_samples.
      #
      # We can now get whatever information we need from the control and signal samples

      # <codecell>

      cut = mctcut
      sel = get_samples( data, cut)
      selzz = get_samples( mczz, cut)
      selz = get_samples( mcz, cut)

      sel['sig_mct_low'+flavor] = sel['sig_mct_low'+flavor] & (data.mctperp < cut)
      sel['sig_mct_low_of'] = sel['sig_mct_low_of'] & (data.mctperp < cut)
      sel['wjets_mct_low'+flavor] = sel['wjets_mct_low'+flavor] & (data.mctperp < cut)
      selzz['sig_mct_low'+flavor] = selzz['sig_mct_low'+flavor] & (mczz.mctperp < cut)
      selz['sig_mct_low'+flavor] = selz['sig_mct_low'+flavor] & (mcz.mctperp < cut)

      # <markdowncell>

      # Get the true number of events from each mc type in the low $M_{\text{CT}}$ control region with b-jet veto and isolation cuts applied

      # <codecell>

      # <markdowncell>

      # Background Normalization
      # ========================
      #
      # Given a shape for a background, all we need to make a prediction for the number of events in the high signal region is some sort of normalization. We get the normalization by fitting for the number of events from each background in the low $M_{CT}$ control region, but which pass all other cuts.
      #
      # Control region:
      #
      # *    low $M_{CT}$: Apply b-tag veto and lepton isolation cuts. $5 < M_{CT} < 100$ GeV
      #
      # The extended PDF for this region is given by adding up the shapes from the control region. Here $S_x(M_{CT})$ indicates the shape as a function of $M_{CT}$ as obtained in control region x
      #
      # $$ P(M_{CT}) = N_{\text{top}}S_{\text{top}}(M_{CT}) + N_{\text{wjets}}S_{\text{wjets}}(M_{CT}) + N_{\text{WW}}S_{\text{WW}}(M_{CT}) $$
      #
      # We use an extended likelihood fit to find $N_{\text{top}}$, $N_{\text{wjets}}$, and $N_{\text{WW}}$.

      # <markdowncell>

      # Load the data from the low $M_{CT}$ control region into a RooDataSet

      # <codecell>

      mct = r.RooRealVar("mct", "mct", 10, mctcut)
      w = r.RooRealVar("w", "w", 0., 10.)
      ds = create_roodataset( data[sel['sig_mct_low'+flavor]].mctperp, data[sel['sig_mct_low'+flavor]].weight, mct, w, title="data")

      # <markdowncell>

      # Build the background shapes using events in the control regions

      # <codecell>

      # first make a TH1
      topshapehistTH1 = create_TH1( data[sel['sig_mct_low_of']].mctperp, data[sel['sig_mct_low_of']].weight, "top", range=(10,mctcut), nbins=int((mctcut-10)/10))
      wjetsshapehistTH1 = create_TH1( data[sel['wjets_mct_low'+flavor]].mctperp, data[sel['wjets_mct_low'+flavor]].weight, "wjets", range=(10,mctcut), nbins=int((mctcut-10)/10))
      vvshapehistTH1 = create_TH1( mczz[selzz['sig_mct_low'+flavor]].mctperp, mczz[selzz['sig_mct_low'+flavor]].weight, "vv", range=(10,mctcut), nbins=int((mctcut-10)/10))
      zshapehistTH1 = create_TH1( mcz[selz['sig_mct_low'+flavor]].mctperp, mcz[selz['sig_mct_low'+flavor]].weight, "Z", range=(10,mctcut), nbins=int((mctcut-10)/10))
      # Then convert to a RooDataHist
      topshapehist = r.RooDataHist("topHist", "topHist", r.RooArgList(mct), topshapehistTH1)
      wjetsshapehist = r.RooDataHist("wjetsHist", "wjetsHist", r.RooArgList(mct), wjetsshapehistTH1)
      vvshapehist = r.RooDataHist("vvHst", "vvHist", r.RooArgList(mct), vvshapehistTH1)
      zshapehist = r.RooDataHist("zHist", "zHist", r.RooArgList(mct), zshapehistTH1)
      # and create a PDF from the RooDataHist
      topshape = r.RooHistPdf( "top", "top", r.RooArgSet(mct), topshapehist)
      wjetsshape = r.RooHistPdf( "wjets", "wjets", r.RooArgSet(mct), wjetsshapehist)
      vvshape = r.RooHistPdf( "vv", "vv", r.RooArgSet(mct), vvshapehist)
      zshape = r.RooHistPdf("z", "z", r.RooArgSet(mct), zshapehist)

      # <codecell>

      # build the extended PDFs
      ntop = r.RooRealVar("ntop", "ntop", 514, 0., 10000.)
      topext = r.RooExtendPdf( "topext", "topext", topshape, ntop)
      nwjets = r.RooRealVar("nwjets", "nwjets", 300, 0., 10000.)
      wjetsext = r.RooExtendPdf( "wjetsext", "wjetsext", wjetsshape, nwjets)
      nvv = r.RooRealVar("nvv", "nvv", 600, 0., 10000.)
      vvext = r.RooExtendPdf( "vvext", "vv", vvshape, nvv)
      nz = r.RooRealVar("nz", "nz", 300, 0, 10000)
      zext = r.RooExtendPdf("zext", "z", zshape, nz)

      # <markdowncell>

      # Add up the shapes to create a composite model

      # <codecell>

      model = r.RooAddPdf("model", "model", r.RooArgList(topext, wjetsext, vvext, zext))

      # <markdowncell>

      # Using b-tagging efficiency to constrain $N_\mathrm{top}$
      # ========================================================
      # Because the top and vv shapes are very similar in the low $M_{CT}$ region, the fit tends to have very large uncertainties on these numbers. We can improve this by providing an indpendent estimate of $N_{\text{top}}$ and putting it in the fit as a Gaussian constraint.
      #
      # First we estimate the b-tagging efficiency. Assume that any event that has at least 1 b-tag actually has 2 b-tags. Let $N_x$ be the number of events with exactly $x$ b-tags and $\epsilon$ be the b-tagging efficiency. Then
      #
      # $$ N_2 = N_{\text{bb}}\epsilon^2 $$
      # $$ N_1 = 2N_{\text{bb}}\epsilon(1-\epsilon) $$
      # $$ N_0 = N_{\text{bb}}(1-\epsilon)^2 $$
      #
      # Then we get
      #
      # $$ \epsilon = \frac{2N_2}{N_1+2N_2} $$
      #
      # with uncertainty
      #
      # $$ \sigma_\epsilon = \epsilon\sqrt{\frac{1}{2N_2}-\frac{1}{N_1+2N_2}} $$

      # <codecell>

      # n_1tag_tt = sum(data[sel['1tag_mct_low']].weight) - sum(mctw[seltw['1tag_mct_low']].weight)
      n_1tag = sum(data[sel['1tag_mct_low'+flavor]].weight)
      n_2tag = sum(data[sel['2tag_mct_low'+flavor]].weight)
      # eff = 2.*n_2tag/(n_1tag_tt+2*n_2tag)
      eff = 2.*n_2tag/(n_1tag+2*n_2tag)
      seff = eff*sqrt(1./2/n_2tag - 1./(n_1tag+2*n_2tag))
      print "Eff: {} +- {}".format(eff, seff)
      corr_eff = 1.057*eff

      # seff = eff*sqrt(1./2/n_2tag - 1./(n_1tag_tt+2*n_2tag))

      # <markdowncell>

      # We can estimate the number of bb+X events with 0 b-tags, which is the number of tt/tW events in our low $M_{CT}$ control region
      #
      # $$N_0 = \frac{N_1}{2}\frac{1-\epsilon}{\epsilon} $$
      #
      # with uncertainty
      #
      # $$\sigma_{N_0} = N_0\sqrt{\frac{2}{N_1}+4\frac{\sigma_\epsilon^2}{\epsilon^2}} $$

      # <codecell>

      ntop_ctrl_high_val = sum(data[sel['top_mct_high'+flavor]].weight)
      ntop_ctrl_high_err_val = scipy.stats.poisson.std(ntop_ctrl_high_val)
      ntop_ctrl_high_var = r.RooRealVar("ntop_ctrl_high", "ntop_ctrl_high", ntop_ctrl_high_val, 0, 10000)
      ntop_ctrl_high_obs = r.RooRealVar("ntop_ctrl_high_obs", "ntop_ctrl_high_obs", ntop_ctrl_high_val)
      ntop_ctrl_high_err = r.RooRealVar("ntop_ctrl_high_err", "ntop_ctrl_high_err", ntop_ctrl_high_err_val)
      ntop_ctrl_high_constraint = r.RooGaussian("ntop_ctrl_high_constraint", "ntop_ctrl_high_constraint", ntop_ctrl_high_var,
                                                 ntop_ctrl_high_obs, ntop_ctrl_high_err)
      ntop_ctrl_low_val = sum(data[sel['top_mct_low'+flavor]].weight)
      ntop_ctrl_low_err_val = scipy.stats.poisson.std(ntop_ctrl_low_val)
      ntop_ctrl_low_var = r.RooRealVar("ntop_ctrl_low", "ntop_ctrl_low", ntop_ctrl_low_val, 0, 10000)
      ntop_ctrl_low_obs = r.RooRealVar("ntop_ctrl_low_obs", "ntop_ctrl_low_obs", ntop_ctrl_low_val)
      ntop_ctrl_low_err = r.RooRealVar("ntop_ctrl_low_err", "ntop_ctrl_low_err", ntop_ctrl_low_err_val)
      ntop_ctrl_low_constraint = r.RooGaussian("ntop_ctrl_low_constraint", "ntop_ctrl_low_constraint", ntop_ctrl_low_var,
                                                 ntop_ctrl_low_obs, ntop_ctrl_low_err)
      nz_ctrl_high_val = sum(mcz[selz['sig_mct_high'+flavor]].weight)
      nz_ctrl_high_err_val = scipy.stats.poisson.std(nz_ctrl_high_val) if nz_ctrl_high_val > 0 else 1
      nz_ctrl_high_var = r.RooRealVar("nz_ctrl_high", "nz_ctrl_high", nz_ctrl_high_val, 0, 10000)
      nz_ctrl_high_obs = r.RooRealVar("nz_ctrl_high_obs", "nz_ctrl_high_obs", nz_ctrl_high_val)
      nz_ctrl_high_err = r.RooRealVar("nz_ctrl_high_err", "nz_ctrl_high_err", nz_ctrl_high_err_val)
      nz_ctrl_high_constraint = r.RooGaussian("nz_ctrl_high_constraint", "nz_ctrl_high_constraint", nz_ctrl_high_var,
                                                 nz_ctrl_high_obs, nz_ctrl_high_err)
      nz_ctrl_low_val = sum(mcz[selz['sig_mct_low'+flavor]].weight)
      nz_ctrl_low_err_val = scipy.stats.poisson.std(nz_ctrl_low_val)
      nz_ctrl_low_var = r.RooRealVar("nz_ctrl_low", "nz_ctrl_low", nz_ctrl_low_val, 0, 10000)
      nz_ctrl_low_obs = r.RooRealVar("nz_ctrl_low_obs", "nz_ctrl_low_obs", nz_ctrl_low_val)
      nz_ctrl_low_err = r.RooRealVar("nz_ctrl_low_err", "nz_ctrl_low_err", nz_ctrl_low_err_val)
      nz_ctrl_low_constraint = r.RooGaussian("nz_ctrl_low_constraint", "nz_ctrl_low_constraint", nz_ctrl_low_var,
                                                 nz_ctrl_low_obs, nz_ctrl_low_err)
      nwjets_ctrl_high_val = sum(data[sel['wjets_mct_high'+flavor]].weight)
      if nwjets_ctrl_high_val == 0:
            nwjets_ctrl_high_err_val = 0.1
      else:
            nwjets_ctrl_high_err_val = scipy.stats.poisson.std(nwjets_ctrl_high_val)
      nwjets_ctrl_high_var = r.RooRealVar("nwjets_ctrl_high", "nwjets_ctrl_high", nwjets_ctrl_high_val, -10000, 10000)
      nwjets_ctrl_high_obs = r.RooRealVar("nwjets_ctrl_high_obs", "nwjets_ctrl_high_obs", nwjets_ctrl_high_val)
      nwjets_ctrl_high_err = r.RooRealVar("nwjets_ctrl_high_err", "nwjets_ctrl_high_err", nwjets_ctrl_high_err_val)
      nwjets_ctrl_high_constraint = r.RooGaussian("nwjets_ctrl_high_constraint", "nwjets_ctrl_high_constraint", nwjets_ctrl_high_var,
                                                 nwjets_ctrl_high_obs, nwjets_ctrl_high_err)
      nwjets_ctrl_low_val = sum(data[sel['wjets_mct_low'+flavor]].weight)
      nwjets_ctrl_low_err_val = scipy.stats.poisson.std(nwjets_ctrl_low_val)
      nwjets_ctrl_low_var = r.RooRealVar("nwjets_ctrl_low", "nwjets_ctrl_low", nwjets_ctrl_low_val, 0, 10000)
      nwjets_ctrl_low_obs = r.RooRealVar("nwjets_ctrl_low_obs", "nwjets_ctrl_low_obs", nwjets_ctrl_low_val)
      nwjets_ctrl_low_err = r.RooRealVar("nwjets_ctrl_low_err", "nwjets_ctrl_low_err", nwjets_ctrl_low_err_val)
      nwjets_ctrl_low_constraint = r.RooGaussian("nwjets_ctrl_low_constraint", "nwjets_ctrl_low_constraint", nwjets_ctrl_low_var,
                                                 nwjets_ctrl_low_obs, nwjets_ctrl_low_err)
      nvv_ctrl_high_val = sum(mczz[selzz['sig_mct_high'+flavor]].weight)
      nvv_ctrl_high_err_val = scipy.stats.poisson.std(nvv_ctrl_high_val)
      nvv_ctrl_high_var = r.RooRealVar("nvv_ctrl_high", "nvv_ctrl_high", nvv_ctrl_high_val, 0, 10000)
      nvv_ctrl_high_obs = r.RooRealVar("nvv_ctrl_high_obs", "nvv_ctrl_high_obs", nvv_ctrl_high_val)
      nvv_ctrl_high_err = r.RooRealVar("nvv_ctrl_high_err", "nvv_ctrl_high_err", nvv_ctrl_high_err_val)
      nvv_ctrl_high_constraint = r.RooGaussian("nvv_ctrl_high_constraint", "nvv_ctrl_high_constraint", nvv_ctrl_high_var,
                                                 nvv_ctrl_high_obs, nvv_ctrl_high_err)
      nvv_ctrl_low_val = sum(mczz[selzz['sig_mct_low'+flavor]].weight)
      nvv_ctrl_low_err_val = scipy.stats.poisson.std(nvv_ctrl_low_val)
      nvv_ctrl_low_var = r.RooRealVar("nvv_ctrl_low", "nvv_ctrl_low", nvv_ctrl_low_val, 0, 10000)
      nvv_ctrl_low_obs = r.RooRealVar("nvv_ctrl_low_obs", "nvv_ctrl_low_obs", nvv_ctrl_low_val)
      nvv_ctrl_low_err = r.RooRealVar("nvv_ctrl_low_err", "nvv_ctrl_low_err", nvv_ctrl_low_err_val)
      nvv_ctrl_low_constraint = r.RooGaussian("nvv_ctrl_low_constraint", "nvv_ctrl_low_constraint", nvv_ctrl_low_var,
                                                 nvv_ctrl_low_obs, nvv_ctrl_low_err)



      # <codecell>

      all_pdfs = r.RooArgList(ntop_ctrl_low_constraint, ntop_ctrl_high_constraint, nz_ctrl_low_constraint, nz_ctrl_high_constraint,
                                     nwjets_ctrl_low_constraint, nwjets_ctrl_high_constraint, nvv_ctrl_low_constraint, nvv_ctrl_high_constraint)
      all_pdfs.add(r.RooArgList(model))
      constraint_vars = r.RooArgSet(ntop_ctrl_low_var, ntop_ctrl_high_var, nz_ctrl_low_var, nz_ctrl_high_var,
                                    nwjets_ctrl_low_var, nwjets_ctrl_high_var, nvv_ctrl_low_var, nvv_ctrl_high_var)

      # <codecell>

      model_c = r.RooProdPdf("model_c", "model_c", all_pdfs)

      # <markdowncell>

      # Do the fit

      # <codecell>

      # turn off logging from the module used to do the multicore pdf evaluation
      # r.RooMsgService.instance().getStream(1).removeTopic(r.RooFit.Minimization)
      # r.RooMsgService.instance().getStream(1).removeTopic(r.RooFit.Eval)

      results = model_c.fitTo( ds,
                               r.RooFit.Constrain(constraint_vars), r.RooFit.Save(),
                               r.RooFit.NumCPU(1), r.RooFit.PrintLevel(1), r.RooFit.PrintEvalErrors(-1),
                               r.RooFit.Verbose(False))

      # <markdowncell>

      # Set up all of the calculated variables for the actual background yields. Doing it this way lets RooFit automagically calculate the errors for us, including the full covariance matrix.

      # <codecell>

      ntop_sig_pred = r.RooFormulaVar("ntop_sig_pred", "ntop_sig_pred", "ntop*ntop_ctrl_high/ntop_ctrl_low",
                                      r.RooArgList(ntop, ntop_ctrl_high_var, ntop_ctrl_low_var))
      nz_sig_pred = r.RooFormulaVar("nz_sig_pred", "nz_sig_pred", "nz*nz_ctrl_high/nz_ctrl_low",
                                      r.RooArgList(nz, nz_ctrl_high_var, nz_ctrl_low_var))
      nvv_sig_pred = r.RooFormulaVar("nvv_sig_pred", "nvv_sig_pred", "nvv*nvv_ctrl_high/nvv_ctrl_low",
                                      r.RooArgList(nvv, nvv_ctrl_high_var, nvv_ctrl_low_var))
      nwjets_sig_pred = r.RooFormulaVar("nwjets_sig_pred", "nwjets_sig_pred", "nwjets*nwjets_ctrl_high/nwjets_ctrl_low",
                                      r.RooArgList(nwjets, nwjets_ctrl_high_var, nwjets_ctrl_low_var))
      ntot_sig_pred = r.RooFormulaVar("ntot_sig_pred", "ntot_sig_pred", "ntop_sig_pred+nz_sig_pred+nvv_sig_pred+nwjets_sig_pred",
                                      r.RooArgList(ntop_sig_pred, nz_sig_pred, nvv_sig_pred, nwjets_sig_pred))


      # <markdowncell>

      # Plot

      # <codecell>

      if plot :

            r.gROOT.SetStyle("Plain")

            f = mct.frame()
            ds.plotOn(f, r.RooFit.Binning(19, 5., 100.))
            model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected))
            model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected), r.RooFit.Components("topext"),\
                r.RooFit.LineColor(2), r.RooFit.LineStyle(r.kDashed), r.RooFit.FillColor(2))
            model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected), r.RooFit.Components("wjetsext"),\
                r.RooFit.LineColor(r.kGreen-1), r.RooFit.LineStyle(r.kDashed))
            model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected), r.RooFit.Components("vvext"),\
                r.RooFit.LineColor(r.kRed+1), r.RooFit.LineStyle(r.kDashed))
            model.plotOn(f, r.RooFit.Normalization( 1.0, r.RooAbsReal.RelativeExpected), r.RooFit.Components("zext"),\
                r.RooFit.LineColor(r.kBlue+1), r.RooFit.LineStyle(r.kDashed))
            f.Draw()
            l1 = r.TLegend(0.1,0.7,0.4,0.9)
            l1.AddEntry("data","Data", "lep")
            e = l1.AddEntry("topext","t\\bar{t} + tW","l")
            e.SetLineColor(2)
            e.SetLineStyle(r.kDashed)
            e.SetLineWidth(3)
            e = l1.AddEntry("wjetsext","W+Jets","l")
            e.SetLineColor(r.kGreen-1)
            e.SetLineStyle(r.kDashed)
            e.SetLineWidth(3)
            e = l1.AddEntry("vvext","vv","l")
            e.SetLineColor(r.kRed+1)
            e.SetLineStyle(r.kDashed)
            e.SetLineWidth(3)
            e=l1.AddEntry("zext","Z/#gamma^{*}","l")
            e.SetLineColor(r.kBlue+1)
            e.SetLineStyle(r.kDashed)
            e.SetLineWidth(3)
            l1.SetFillColor(r.kWhite)
            l1.Draw()

      # <markdowncell>

      # Get the values out of the fit and compare to MC values

      # <codecell>

      # ntop_low_mct = ntop.getVal()
      # ntop_low_mct_err = ntop.getError()
      # nwjets_low_mct = nwjets.getVal()
      # nwjets_low_mct_err = nwjets.getError()
      # nvv_low_mct = nvv.getVal()
      # nvv_low_mct_err = nvv.getError()
      # nz_low_mct = nz.getVal()
      # nz_low_mct_err = nz.getError()
      # nww_low_mct = nww.getPropagatedError(results)

      pred_high_dict = {}
      pred_high_dict['Total'] = (ntot_sig_pred.getVal(), ntot_sig_pred.getPropagatedError(results))
      pred_high_dict['Top'] = (ntop_sig_pred.getVal(), ntop_sig_pred.getPropagatedError(results))
      pred_high_dict['VV'] = (nvv_sig_pred.getVal(), nvv_sig_pred.getPropagatedError(results))
      pred_high_dict['DY'] = (nz_sig_pred.getVal(), nz_sig_pred.getPropagatedError(results))
      pred_high_dict['W'] = (nwjets_sig_pred.getVal(), nwjets_sig_pred.getPropagatedError(results))

      pred_low_dict = {}
      pred_low_dict['Top'] = (ntop.getVal(), ntop.getError())
      pred_low_dict['VV'] = (nvv.getVal(), nvv.getError())
      pred_low_dict['DY'] = (nz.getVal(), nz.getError())
      pred_low_dict['W'] = (nwjets.getVal(), nwjets.getError())

      result = {'low':pred_low_dict, 'high':pred_high_dict}
      fout = open("results"+flavor+".json", 'w')
      json.dump(result, fout)
      fout.close()

      return result


