from general import *
from math import sqrt
import ROOT
import pylab
from CMSPyLibs.rootplot import *
from data_handling import *
from filenames import *

"""
Functions to plot variables. They can be plotted using either ROOT or matplotlib. Switch using
set_plot_ROOT() or set_plot_mpl().
"""

ROOT.gStyle.SetPalette(1)

def mpl_hist_func(*args, **kwargs) :
    """
    Wrapper function for pylab.hist that takes the same arguments as plot_hist_ROOT
    """
    name = ""
    if 'name' in kwargs.keys() :
        name = kwargs['name']
        del kwargs['name']
    pylab.hist(*args, **kwargs)
    pylab.xlabel(name)

hist_func = mpl_hist_func
graph_func = pylab.scatter

def set_plot_ROOT() :
    global hist_func, graph_func
    hist_func = plot_hist_ROOT
    graph_func = plot_graph_ROOT

def set_plot_mpl() :
    global hist_func, graph_func
    hist_func = mpl_hist_func
    graph_func = pylab.scatter

def plot_stacked_mct_hist(bins=20, histrange=(0,300)) :
    """Plot a nice stacked histogram of the mct distributions for the signal and background"""
    lm6_mct   = [x["mct"] for x in load_data_iterator_cut(lm6_hd5, """(relIso2 < 0.17) & (mct > 0.5)""")]
    tt_mct    = [x["mct"] for x in load_data_iterator_cut(tt_hd5, """(relIso2 < 0.17) & (mct > 0.5)""")]
    wjets_mct = [x["mct"] for x in load_data_iterator_cut(wjets_hd5, """(relIso2 < 0.17) & (mct > 0.5)""")]

    lm6_hist = ROOT.TH1D("LM6", "LM6", bins, histrange[0], histrange[1])
    tt_hist = ROOT.TH1D("ttJets", "ttJets", bins, histrange[0], histrange[1])
    wjets_hist = ROOT.TH1D("WJets", "WJets", bins, histrange[0], histrange[1])

    for data,hist, weight in zip([lm6_mct, tt_mct, wjets_mct],[lm6_hist, tt_hist, wjets_hist], [lm6_w, tt_w, wjets_w]) :
        for d in data :
            hist.Fill(d, weight)

    lm6_hist.SetFillColor(ROOT.kBlue)
    tt_hist.SetFillColor(ROOT.kRed+1)
    wjets_hist.SetFillColor(ROOT.kGreen-1)
    lm6_hist.SetLineColor(ROOT.kBlue)
    tt_hist.SetLineColor(ROOT.kRed+1)
    wjets_hist.SetLineColor(ROOT.kGreen-1)

    stack = ROOT.THStack("stack", "MCTPerp")
    stack.Add(wjets_hist)
    stack.Add(tt_hist)
    stack.Add(lm6_hist)

    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptLogy()

    stack.Draw()

    l = ROOT.TLegend(0.7,0.7,0.9,0.9)
    l.AddEntry("LM6", "LM6", "f")
    l.AddEntry("ttJets", "TTJets", "f")
    l.AddEntry("WJets", "WJets", "f")
    l.SetFillColor(ROOT.kWhite)
    l.Draw()
    raw_input("...")

def plot_data_hist( hdf_file, data_item, bins=30, range=(0,300)) :
    """
    plot a distribution from an hdf file. data_item should be a string containing the
    name of the field
    """
    data = pylab.array([x[data_item] for x in load_data_iterator_iso(hdf_file)])
    h = hist_func(data, bins, range=range, name=data_item)
    return h

def plot_ptmax_hist( hdf_file, bins=30, range=(0,300)) :
    """
    plot a distribution from an hdf file. data_item should be a string containing the
    name of the field
    """
    pts = [(x["pt1"],x["pt2"]) for x in load_data_iterator_iso(hdf_file)]
    ptmax = map(max, pts)
    h = hist_func(ptmax, bins, range=range, name="ptLeadingLepton")
    return h

def plot_ptmin_hist( hdf_file, bins=30, range=(0,300)) :
    """
    plot a distribution from an hdf file. data_item should be a string containing the
    name of the field
    """
    pts = [(x["pt1"],x["pt2"]) for x in load_data_iterator_iso(hdf_file)]
    ptmax = map(min, pts)
    h = hist_func(ptmax, bins, range=range, name="ptLepton2")
    return h

def plot_ptmax_vs_ptmin( hdf_file, bins=30, range=(0,300)) :
    """
    plot a distribution from an hdf file. data_item should be a string containing the
    name of the field
    """
    pts = [(x["pt1"],x["pt2"]) for x in load_data_iterator_iso(hdf_file)]
    ptmax = map(max, pts)
    ptmin = map(min, pts)
    h = plot_graph_ROOT(ptmax, ptmin)
    return h

def plot_sqrtptprod_hist( hdf_file, bins=30, range=(0,300)) :
    """
    plot a distribution from an hdf file. data_item should be a string containing the
    name of the field
    """
    ptprod = [sqrt(x["pt1"]*x["pt2"]) for x in load_data_iterator_iso(hdf_file)]
    h = hist_func(ptprod, bins, range=range, name="sqrtpt1pt2")
    return h

def plot_phiUp_hist( hdf_file, bins=30, range=(0,3.14)) :
    """
    plot a distribution from an hdf file. data_item should be a string containing the
    name of the field
    """
    phiUp = [x["phiUp"] for x in load_data_iterator_iso(hdf_file)]
    h = hist_func(phiUp, bins, range=range, name="phiUp")
    return h

def plot_phiUp_vs_delta( hdf_file) :
    """
    plot a distribution from an hdf file. data_item should be a string containing the
    name of the field
    """
    data = [(x["phiUp"],x["delta"]) for x in load_data_iterator_iso(hdf_file)]
    phiUp, delta = zip(*data)
    phiUp2 = []
    for phi in phiUp :
        if phi < pylab.pi :
            phiUp2.append(phi)
        else :
            phiUp2.append(2*pylab.pi - phi)
    h = plot_hist2d_ROOT(phiUp2, delta, xbins=5, ybins=5, xrange=(0,3.14), yrange=(0, 3.14/2))
    return h

def plot_ptprod_vs_phiUp( hdf_file) :
    """
    plot a distribution from an hdf file. data_item should be a string containing the
    name of the field
    """
    data = [(sqrt(x["pt1"]*x["pt2"]),x["phiUp"]) for x in load_data_iterator_iso(hdf_file)]
    ptprod, phiUp = zip(*data)
    phiUp2 = []
    for phi in phiUp :
        if phi < pylab.pi :
            phiUp2.append(phi)
        else :
            phiUp2.append(2*pylab.pi - phi)
    # h = plot_hist2d_ROOT(ptprod, phiUp2, xbins=5, ybins=5, xrange=(0,150.), yrange=(0, 3.14))
    h = graph_func(ptprod, phiUp2)
    return h

def plot_ptprod_vs_delta( hdf_file, bins=30, range=(0,3.14)) :
    """
    plot a distribution from an hdf file. data_item should be a string containing the
    name of the field
    """
    data = [(sqrt(x["pt1"]*x["pt2"]),x["delta"]) for x in load_data_iterator_iso(hdf_file)]
    ptprod, delta = zip(*data)
    h = graph_func(ptprod, delta)
    return h








