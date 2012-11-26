"""Check the closure of the top veto efficiency estimate"""

n_1tag = sum(mc[smc['1tag_mct_low']].weight)
n_2tag = sum(mc[smc['2tag_mct_low']].weight)
eff = 2.*n_2tag/(n_1tag+2*n_2tag)
seff = eff*sqrt(1./2/n_2tag - 1./(n_1tag+2*n_2tag))
      # $$\sigma_{N_0} = N_0\sqrt{\frac{2}{N_1}+4\frac{\sigma_\epsilon^2}{\epsilon^2}} $$
print "Eff: {} +- {}".format(eff, seff)

ntop_pred = n_1tag/2.*(1-eff)/eff
ntop_pred_err = ntop_pred*sqrt(2./n_1tag+4*seff**2/eff**2)
ntop_true = sum(mc[smc['sig_mct_low']&(mc.mc_cat=='top')].weight)
print "True: {}\tPredicted: {} +- {}".format(ntop_true, ntop_pred, ntop_pred_err)

