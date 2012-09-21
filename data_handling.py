import pandas
import logging
import numpy
import ROOT
import itertools
import functools
import numpy.random
import CMSPyLibs.events
import CMSPyLibs.general_calc as general_calc
import CMSPyLibs.event_counter as event_counter
from CMSPyLibs.cmsutilities import get_TLorentzVector, get_5_X_isolation, angle_0_2pi, get_upstream_phi_res

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
    getter = CMSPyLibs.events.CMSDileptonEventGetter(input_files)
    getter.set_electron_collection("nonIsolatedPatElectrons")
    getter.set_muon_collection("nonIsolatedPatMuons")
    getter.set_jet_collection("selectedPatJets")
    getter.vertex_collection = "offlinePrimaryVertices"
    getter.set_met_collection("patMETs")
    getter.do_PU = (mctype != 'data')
    getter.do_SMS = ('sms' in mctype)

    count = event_counter.EventCounter()
    calc = general_calc.VarCalculator()

    datadicts = []

    eventids = set()

    for event in getter.events() :
        count.Increment()

        if event.eventID in eventids :
            print "Duplicate event, skipping..."
        else :
            eventids.add(event.eventID)
        event.jet_btag = "combinedSecondaryVertexBJetTags"
        event.jet_btag_cut = 0.679
        datarow = {}
        datarow["mctype"] = mctype
        datarow['run'] = event.eventID.run_number
        datarow['lumi'] = event.eventID.luminosity_block
        datarow['event'] = event.eventID.event_number

        if mctype == 'data':
            datarow['DoubleEle_Trigger'] = event.passes_HLT("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")
            datarow['DoubleMu_Trigger'] = event.passes_HLT("HLT_Mu17_Mu8_v*")
            datarow['EMu_Trigger'] = event.passes_HLT("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")\
                                    or event.passes_HLT("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")

        datarow['nvertices'] =len(event.get_vertices())

        if 'num_pu_vertices' in event.metadata.keys() :
            datarow['nPuVertices'] = event.metadata['num_pu_vertices']

        if 'modelParams' in event.metadata.keys() :
            datarow['mass1'] = event.metadata['modelParams'][0]
            datarow['mass2'] = event.metadata['modelParams'][1]

        # lepton info
        get_isolation = functools.partial(get_5_X_isolation, rho = event.metadata['rho'])
        leptons = event.get_leptons()
        leptons.sort(key = get_isolation)
        datarow["relIso1"] = get_isolation( leptons[0] )
        datarow["relIso2"] = get_isolation( leptons[1] )
        datarow["pt1"] = leptons[0].pt()
        datarow["pt2"] = leptons[1].pt()
        datarow["pdg1"] = leptons[0].pdgId()
        datarow["pdg2"] = leptons[1].pdgId()
        datarow["phi1"] = leptons[0].phi()
        datarow["phi2"] = leptons[1].phi()
        datarow["eta1"] = leptons[0].eta()
        datarow["eta2"] = leptons[1].eta()
        try :
            datarow['parentPdg1'] = get_parent(leptons[0].genParticle()).pdgId()
            datarow['parentPdg2'] = get_parent(leptons[1].genParticle()).pdgId()
        except ReferenceError :
            pass
        # is there a third isolated lepton?
        datarow['ThirdLepton'] = False
        if len(leptons) > 2:
            if get_isolation(leptons[2]) < 0.15:
                datarow['ThirdLepton'] = True

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
        if twodelta > numpy.pi : twodelta = 2*numpy.pi-twodelta
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
        # if event.get_met().pt() < 30 :
        #     logging.debug("Skipping event with MET "+str(event.get_met().pt()))
        #     continue
        datarow['metPhi'] = event.get_met().phi()
        datarow['weight'] = weight

        # invariant mass
        datarow["mll"] = (get_TLorentzVector(leptons[0])+get_TLorentzVector(leptons[1])).M()

        logging.debug("Adding event")
        datadicts.append(datarow)

    logging.info("Saving to "+output_file)
    store = pandas.HDFStore(output_file)
    df = pandas.DataFrame(datadicts)
    if mctype in store.keys() :
	logging.debug( mctype+" already exists. Appending.")
        store[mctype] = store[mctype].append(df, ignore_index=True)
    else :
	logging.debug("Creating new DataFrame")
        store.put(mctype, df)
    store.flush()
    store.close()
