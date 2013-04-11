from collections import defaultdict
import ROOT

ROOT.gROOT.ProcessLineSync('.L IntegralError.C+')

dummy_stat = ROOT.RooStats.NumEventsTestStat

def get_results(ws, res):

    results = defaultdict(lambda: defaultdict(dict))

    # ROOT.gROOT.ProcessLineSync(".L IntegralError.C+")
    obs_sf = ws.obj("obs_x_sf")
    obs_of = ws.obj("obs_x_of")

    tf, cov = get_tf1("sf_model", ws, res, obs_sf)
    results['sf']['sum'] = get_integrals(tf, cov, 1.)

    tf, cov = get_tf1("top_sf_sf_overallSyst_x_StatUncert", ws, res, obs_sf)
    results['sf']['top'] = get_integrals(tf, cov)

    tf, cov = get_tf1("vv_sf_sf_overallSyst_x_StatUncert_x_sf_ww_syst_sf_ShapeSys", ws, res, obs_sf)
    results['sf']['vv'] = get_integrals(tf, cov)

    tf, cov = get_tf1("z_sf_sf_overallSyst_x_StatUncert_x_sf_z_syst_sf_ShapeSys", ws, res, obs_sf)
    results['sf']['z'] = get_integrals(tf, cov)

    tf, cov = get_tf1("wjets_sf_sf_overallSyst_x_StatUncert_x_sf_wjets_syst_sf_ShapeSys", ws, res, obs_sf)
    results['sf']['fake'] = get_integrals(tf, cov)

    tf, cov = get_tf1("of_model", ws, res, obs_of)
    results['of']['sum'] = get_integrals(tf, cov, 1.)

    tf, cov = get_tf1("top_of_of_overallSyst_x_StatUncert", ws, res, obs_of)
    results['of']['top'] = get_integrals(tf, cov)

    tf, cov = get_tf1("vv_of_of_overallSyst_x_StatUncert_x_of_ww_syst_of_ShapeSys", ws, res, obs_of)
    results['of']['vv'] = get_integrals(tf, cov)

    tf, cov = get_tf1("z_of_of_overallSyst_x_StatUncert_x_of_z_syst_of_ShapeSys", ws, res, obs_of)
    results['of']['z'] = get_integrals(tf, cov)

    tf, cov = get_tf1("wjets_of_of_overallSyst_x_StatUncert_x_of_wjets_syst_of_ShapeSys", ws, res, obs_of)
    results['of']['fake'] = get_integrals(tf, cov)  

    return results

def get_integrals(tf, cov, multiplier=1./10):
    low, high, up = 10, 120, 300

    low_i = tf.Integral(low, high)*multiplier
    low_err = ROOT.get_integral_error(tf, low, high, cov)*multiplier
    high_i = tf.Integral(high, up)*multiplier
    high_err = ROOT.get_integral_error(tf, high, up, cov)*multiplier

    return {'low':(low_i, low_err),'high':(high_i, high_err)}

def get_tf1(name, ws, res, obs):
    shape = ws.obj(name)
    params = shape.getParameters(ROOT.RooArgSet(obs))
    ROOT.RooStats.RemoveConstantParameters(params)
    cov = res.reducedCovarianceMatrix(ROOT.RooArgList(params))

    return shape.asTF(ROOT.RooArgList(obs), ROOT.RooArgList(params)), cov