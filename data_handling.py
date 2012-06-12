import pandas
import pylab
import ROOT
import itertools
import numpy.random
import CMSPyLibs.events.dilepton_event as dilepton_event
import CMSPyLibs.general_calc as general_calc
import CMSPyLibs.event_counter as event_counter
from CMSPyLibs.cmsutilities import get_TLorentzVector, get_PF_isolation, angle_0_2pi, get_upstream_phi_res

"""
Functions for handling fit data in the pytables hdf5 format
Includes creating the files and different ways to iterate over data
in a file.
"""

def get_parent( gen_particle ) :
    """Get the parent particle that comes directly from a vector boson"""
    if abs(gen_particle.mother().pdgId()) in range(22,25) :
        return gen_particle
    else :
        return get_parent(gen_particle.mother())

def save_data_pandas( input_files, output_file, mctype="mc", weight=1.):
    """Load all of the relevant data from the CMSSW files and save it using pandas"""
    getter = dilepton_event.CMSDileptonEventGetter(input_files)
    getter.set_electron_collection("loosePatElectrons")
    getter.set_muon_collection("loosePatMuons")

    count = event_counter.EventCounter()
    calc = general_calc.VarCalculator()

    datadicts = []

    for event in getter.events() :
        event.jet_btag = "combinedSecondaryVertexBJetTags"
        event.jet_btag_cut = 0.244
        datarow = {}
        datarow["mctype"] = mctype
        count.Increment()

        # lepton info
        leptons = event.get_leptons()
        leptons.sort(key = get_PF_isolation)
        datarow["relIso1"] = get_PF_isolation( leptons[0])
        datarow["relIso2"] = get_PF_isolation( leptons[1])
        datarow["pt1"] = leptons[0].pt()
        datarow["pt2"] = leptons[1].pt()
        datarow["pdg1"] = leptons[0].pdgId()
        datarow["pdg2"] = leptons[1].pdgId()
        datarow["phi1"] = leptons[0].phi()
        datarow["phi2"] = leptons[1].phi()
        try :
            datarow['parentPdg1'] = get_parent(leptons[0].genParticle()).pdgId()
            datarow['parentPdg2'] = get_parent(leptons[1].genParticle()).pdgId()
        except ReferenceError :
            pass


        # jets
        datarow['njets'] = len(event.get_jets())
        datarow['nbjets'] = len(event.get_tagged_jets())

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
        datarow["mct"], datarow['mct_type'] = calc.mctBoostCor_210()
        datarow['mctperp'] = calc.mctPerp_210()

        # met
        datarow["metPt"] = event.get_met().pt()
        datarow['metPhi'] = event.get_met().phi()
        datarow['weight'] = weight

        # invariant mass
        datarow["mll"] = (get_TLorentzVector(leptons[0])+get_TLorentzVector(leptons[1])).M()

        datadicts.append(datarow)

    store = pandas.HDFStore(output_file)
    df = pandas.DataFrame(datadicts)
    if mctype in store.keys() :
        store[mctype] = store[mctype].append(df, ignore_index=True)
    else :
        store.put(mctype, df)
    store.flush()
    store.close()
