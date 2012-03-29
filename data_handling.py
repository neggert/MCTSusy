import tables
import pylab
import ROOT
import CMSPyLibs.events.dilepton_event as dilepton_event
import CMSPyLibs.general_calc as general_calc
import CMSPyLibs.event_counter as event_counter
from CMSPyLibs.cmsutilities import get_TLorentzVector, get_PF_isolation, angle_0_2pi, get_upstream_phi_res

"""
Functions for handling fit data in the pytables hdf5 format
Includes creating the files and different ways to iterate over data
in a file.
"""

class FitEvent(tables.IsDescription):
    idnumber = tables.Int64Col()
    pt1      = tables.Float64Col()
    pt2      = tables.Float64Col()
    met      = tables.Float64Col()
    delta    = tables.Float64Col()
    phiUp    = tables.Float64Col()
    sPhiUp   = tables.Float64Col()
    relIso1  = tables.Float64Col()
    relIso2  = tables.Float64Col()
    mct      = tables.Float64Col()
    mctype   = tables.StringCol(8)

def save_data_hdf5( input_files, output_file, mctype="mc" ):
    """Load all of the relevant data from the CMSSW files and save it using pyTables"""
    getter = dilepton_event.CMSDileptonEventGetter(input_files)
    getter.set_electron_collection("loosePatElectrons")
    getter.set_muon_collection("loosePatMuons")

    h5file = tables.openFile(output_file, "w", "Fit Data File")
    group = h5file.createGroup("/", "fit", "Fit information")

    table = h5file.createTable( group, "fitdata", FitEvent, "Fit Data")

    count = event_counter.EventCounter()
    calc = general_calc.VarCalculator()

    try :
        for event in getter.events() :
            datarow = table.row
            datarow["idnumber"] = count.eventNum
            datarow["mctype"] = mctype
            count.Increment()

            # lepton info
            leptons = event.get_leptons()
            leptons.sort(key = get_PF_isolation)
            datarow["relIso1"] = get_PF_isolation( leptons[0])
            datarow["relIso2"] = get_PF_isolation( leptons[1])
            datarow["pt1"] = leptons[0].pt()
            datarow["pt2"] = leptons[1].pt()

            # angles
            phi1 = leptons[0].phi()
            phi2 = leptons[1].phi()
            phiUp = event.upstream_of_leptons().Phi()
            phiUpSigMatrix = event.get_met().getSignificanceMatrix()
            gamma = angle_0_2pi((phi1+phi2)/2) # average of visible angles
            twodelta = abs(phi1-phi2)
            if twodelta > pylab.pi : twodelta = 2*pylab.pi-twodelta
            datarow["delta"] = twodelta/2
            datarow["phiUp"] = angle_0_2pi(phiUp - gamma)
            datarow["sPhiUp"] = get_upstream_phi_res(event.upstream_of_leptons(), phiUpSigMatrix)

            # mct
            p4jet1 = ROOT.TLorentzVector(0., 0., 0., 0.)
            p4jet2 = ROOT.TLorentzVector(0., 0., 0., 0.)
            met = event.get_met()
            calc.setP4s(p4jet1, p4jet2, get_TLorentzVector(leptons[0]), get_TLorentzVector(leptons[1]), get_TLorentzVector(met))
            datarow["mct"] = calc.mctPerp_210()

            datarow["met"] = event.get_met().pt()

            datarow.append()
    finally : 
        table.flush()
        h5file.close()

def load_data_iterator( hd5_file ) :
    """Load the fit data from the pytables file"""
    f = tables.openFile(hd5_file, "r")
    table = f.root.fit.fitdata
    return table.iterrows()

def load_data_iterator_cut( hd5_file, cut ) :
    """Load the fit data from the pytables file with a cut specified by the string cut"""
    f = tables.openFile(hd5_file, "r")
    table = f.root.fit.fitdata
    return table.where(cut)

def load_data_iterator_iso( hd5_file ) :
    """Load the fit data from the pickle file with an isolation cut"""
    f = tables.openFile(hd5_file, "r")
    table = f.root.fit.fitdata
    return table.where("""relIso2 < 0.17""")