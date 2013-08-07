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
from math import sqrt, isnan

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load("libCondFormatsJetMETObjects.so")

from ROOT import std

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

def save_data_pandas( input_files, output_file, mctype="mc", mc_cat="mc_cat", x_eff=1.):
    """Load all of the relevant data from the CMSSW files and save it using pandas"""
    getter = CMSPyLibs.events.CMSEventGetter(input_files)
    getter.set_electron_collection("nonIsolatedPatElectrons")
    getter.set_muon_collection("nonIsolatedPatMuons")
    getter.set_jet_collection("selectedPatJets")
    getter.vertex_collection = "offlinePrimaryVertices"
    getter.set_met_collection("patMETs")
    getter.do_PU = ('data' not in mctype)
    getter.do_genparticles = ('data' not in mctype)

    getter.do_SMS = ('sms' in mctype)

    jec_parameters = ROOT.JetCorrectorParameters("/home/uscms33/MCTSusy/JEC/Summer12_V2_DATA_AK5PF_UncertaintySources.txt", "SubTotalDataMC")
    jec_uncertainty = ROOT.JetCorrectionUncertainty(jec_parameters)

    jec_resid_parameters = ROOT.JetCorrectorParameters("/home/uscms33/MCTSusy/JEC/GR_R_52_V9_L2L3Residual_AK5PF.txt")
    jec_resid_parameters_v = std.vector(ROOT.JetCorrectorParameters)(0)
    jec_resid_parameters_v.push_back(jec_resid_parameters)

    jec_resid = ROOT.FactorizedJetCorrector(jec_resid_parameters_v)

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

        if 'data' in mctype:
            datarow['DoubleEle_Trigger'] = event.passes_HLT("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")
            datarow['DoubleMu_Trigger'] = event.passes_HLT("HLT_Mu17_Mu8_v*")
            datarow['EMu_Trigger'] = event.passes_HLT("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")\
                                    or event.passes_HLT("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")

        datarow['nvertices'] =len(event.get_vertices())

        if 'num_pu_vertices' in event.metadata.keys() :
            datarow['nPuVertices'] = event.metadata['num_pu_vertices']
            datarow['nTruePuVertices'] = event.metadata['num_true_pu_vertices']

        if 'modelParams' in event.metadata.keys() :
            datarow['mass1'] = event.metadata['modelParams']
            #datarow['mass2'] = event.metadata['modelParams'][1]

        # generator level info
        if mctype != 'data':
            neutrinos = [12, 14, 16]
            datarow['gen_neutrinos'] = 0
            for n in event.get_genparticles():
                if abs(n.pdgId()) in neutrinos and n.status()==3:
                    datarow['gen_neutrinos'] += 1

        # # lepton info
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
        try :
            datarow['parentParentPdg1'] = get_parent(leptons[0].genParticle()).mother().pdgId()
            datarow['parentParentPdg2'] = get_parent(leptons[1].genParticle()).mother().pdgId()
        except ReferenceError :
            pass
        # is there a third isolated lepton?
        datarow['ThirdLepton'] = False
        if len(leptons) > 2:
            if get_isolation(leptons[2]) < 0.15:
                datarow['ThirdLepton'] = True
                datarow['relIso3'] = get_isolation(leptons[2])
                datarow['pt3'] = leptons[2].pt()
                datarow['pdg3'] = leptons[2].pdgId()
                datarow['phi3'] = leptons[2].phi()
                datarow['eta3'] = leptons[2].eta()

        # jets
        datarow['njets'] = len(event.get_jets())
        datarow['nbjets'] = len(event.get_tagged_jets())

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
        datarow['x_eff'] = float(x_eff)
        datarow['mc_cat'] = mc_cat

        # invariant mass
        datarow["mll"] = (get_TLorentzVector(leptons[0])+get_TLorentzVector(leptons[1])).M()

        # MET Uncertainty stuff
        if "sms" in mctype:
            scaled_mets = [datarow["metPt"]]
            scaled_mctperps = [datarow["mctperp"]]

            p4met = get_TLorentzVector(met)
            # loop through jets
            p4met_up = p4met
            p4met_down = p4met
            for jet in event.get_jets():
                # get uncorrected jet
                jet_u = get_TLorentzVector(jet.correctedJet(0))

                # add it to the MET vector
                p4met_up += jet_u
                p4met_down += jet_u

                # get its uncertainty
                # want to use corrected kinematics here
                jec_uncertainty.setJetPt(jet.pt())
                jec_uncertainty.setJetEta(jet.eta())
                jec_unc_val = jec_uncertainty.getUncertainty(True)

                # get its residual correction
                # use uncorrected kinematics
                jec_resid.setJetPt(jet_u.Pt())
                jec_resid.setJetEta(jet_u.Eta())
                jec_resid_val = abs(1-jec_resid.getCorrection())

                # add them in quadrature
                unc = sqrt(jec_unc_val**2+jec_resid_val**2)

                # scale the jet by that amount
                jet_scaleup = jet_u*(1+unc)
                jet_scaledown = jet_u*(1-unc)

                # subtract it back from the MET
                p4met_up -= jet_scaleup
                p4met_down -= jet_scaledown

            # calculate MCTPerp with the modified MET
            calc.setP4s(p4jet1, p4jet2, get_TLorentzVector(leptons[0]), get_TLorentzVector(leptons[1]), p4met_up)
            up = calc.mctPerp_210()
            calc.setP4s(p4jet1, p4jet2, get_TLorentzVector(leptons[0]), get_TLorentzVector(leptons[1]), p4met_down)
            down = calc.mctPerp_210()

            scaled_mctperps.extend([up, down])

            scaled_mets.extend([p4met_up.Pt(), p4met_down.Pt()])

            # now move on to unclustered met
            # add up the clustered pt
            clustered_p4 = ROOT.TLorentzVector(0., 0., 0., 0.)
            for jet in event.get_jets():
                # get uncorrected jet
                clustered_p4 += get_TLorentzVector(jet.correctedJet(0))

            for lepton in event.get_leptons():
                clustered_p4 += get_TLorentzVector(lepton)

            p4met_without_clustered = p4met+clustered_p4

            # scale what's left by 10% and subtract the clustered met
            p4met_up = p4met_without_clustered*1.1-clustered_p4
            p4met_down = p4met_without_clustered*0.9-clustered_p4

            # calculate MCTPerp with the modified MET
            calc.setP4s(p4jet1, p4jet2, get_TLorentzVector(leptons[0]), get_TLorentzVector(leptons[1]), p4met_up)
            up = calc.mctPerp_210()
            calc.setP4s(p4jet1, p4jet2, get_TLorentzVector(leptons[0]), get_TLorentzVector(leptons[1]), p4met_down)
            down = calc.mctPerp_210()

            scaled_mctperps.extend([up, down])

            scaled_mets.extend([p4met_up.Pt(), p4met_down.Pt()])

            datarow['mctperp_up'] = max(scaled_mctperps)
            datarow['mctperp_down'] = min(scaled_mctperps)
            datarow['metPt_up'] = max(scaled_mets)
            datarow['metPt_down'] = min(scaled_mets)

        logging.debug("Adding event")
        datadicts.append(datarow)

    logging.info("Saving to "+output_file)
    store = pandas.HDFStore(output_file)
    df = pandas.DataFrame(datadicts)
    if "data" in store.keys() :
	logging.debug( "data already exists. Appending.")
        store['data'] = store["data"].append(df, ignore_index=True)
    else :
	logging.debug("Creating new DataFrame")
        store.put("data", df)
    store.flush()
    store.close()
