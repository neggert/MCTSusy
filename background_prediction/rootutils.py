"""Utilities for creating ROOT and RooFit objects"""

import ROOT as r

def create_roodataset( data, weights, roovar, rooweightvar, title="dataset") :
    """Create a single-variable RooDataSet in roovar"""
    ds = r.RooDataSet( title, title, r.RooArgSet(roovar, rooweightvar), rooweightvar.GetName())
    for d, w in zip(data, weights) :
        roovar.setVal(d)
        ds.add(r.RooArgSet(roovar), w)
    return ds

def create_TH1( data, weights, name="hist") :
    h = r.TH1D(name, name, 19, 5., 100.)
    for d, w in zip(data, weights) :
        h.Fill(d, w)
    return h