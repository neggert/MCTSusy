
def get_samples( data, mctcut=100., real_data=False) :
    """
    Get different control and signal events in sample data.

    Usage:
        # data is of type pandas.DataFrame
        selections = get_samples(data)
        sig_events = data[selections['signal']]

    *data* A pandas.DataFrame instance

    --------------------

    returns: A dict of boolean selectors selecting events that pass
        the selection described by the key. Dict keys and the corresponding
        selections are:

        # isolation-related
        'l1_passes_tight_iso': Lepton 1 passes isolation cuts
        'l2_passes_tight_iso': Lepton 2 passes isolation cuts
        'l1_passes_loose_iso': Lepton 1 passes loose isolation cuts
        'l2_passes_loose_iso': Lepton 2 passes loose isolation cuts
        'l1_passes_ctrl_iso': Lepton 1 falls in iso control region. It passes loose isolation cuts, but not tight
        'l2_passes_ctrl_iso': Lepton 2 falls in iso control region. It passes loose isolation cuts, but not tight
        'isolation_ctrl': 1 leptons passes tight isolation and 1 lepton falls in iso control region
        'isolation_sig': Both leptons pass tight isolation

        # b-tag-related
        '0_bjets' : 0 b-jets
        '1_bjets' : 1 b-jet
        '2_bjets' : 2 b-jets
        'more_bjets' : >2 b-jets
        'bjets_sig' : no b-jets
        'bjets_ctrl': >= 1 b-jet

        # combine iso and b-tag to create control samples
        'wjets_ctrl' : 0 b-jets. 1 control and 1 tight lepton
        '1tag_ctrl' : 1 b-jet. 2 tight leptons
        '2tag_ctrl' : 2 b-jets. 2 tight leptons
        'top_ctrl' : >=1 b-jet. 2 tight leptons
        'iso_bjet_sig' : 0 b-jets. 2 tight leptons

        # mct-related
        'mct_low' : 5 < mct < 100
        'mct_high' : mct > 100

        # divide control samples by mct
        'sig_mct_low': passes iso_bjet_sig and mct_low
        'sig_mct_high' : passes iso_bjet_sig and mct_high
        'wjets_mct_low' : passes wjets_ctrl and mct_low
        'wjets_mct_high' : passes wjets_ctrl and mct_high
        'top_mct_low' : passes top_ctrl and mct_low
        'top_mct_high' : passes top_ctrl and mct_high
        '1tag_mct_low' : passes 1tag_ctrl and mct_low
        '2tag_mct_low' :passes 2tag_ctrl and mct_low


    """
    outdict = {}

    if real_data:
        outdict['ee'] = (abs(data.pdg1) == 11) & (abs(data.pdg2) == 11) & (data.DoubleEle_Trigger)
        outdict['mumu'] = (abs(data.pdg1) == 13) & (abs(data.pdg2) == 13) & (data.DoubleMu_Trigger)
        outdict['emu'] = (abs(data.pdg1) != abs(data.pdg2)) & (data.EMu_Trigger)
    else :
        outdict['ee'] = ((abs(data.pdg1) == 11) & (abs(data.pdg2) == 11))
        outdict['mumu'] = ((abs(data.pdg1) == 13) & (abs(data.pdg2) == 13))
        outdict['emu'] = (abs(data.pdg1) != abs(data.pdg2))
    outdict['opposite_sign'] = (data.pdg1*data.pdg2 < 0)
    outdict['same_sign'] = (data.pdg1*data.pdg2 > 0)  & (~data.ThirdLepton) & (data.mll > 12)
    outdict['opposite_sign_ee'] = outdict['opposite_sign'] & outdict['ee']
    outdict['opposite_sign_mumu'] = outdict['opposite_sign'] & outdict['mumu']
    outdict['opposite_sign_emu'] = outdict['opposite_sign'] & outdict['emu']
    outdict['preselection'] = (outdict['opposite_sign_ee'] | outdict['opposite_sign_mumu'] | outdict['opposite_sign_emu']) & (~data.ThirdLepton) & (data.mll > 12)
    outdict['3leptons'] = outdict['opposite_sign'] & (data.ThirdLepton) & (data.mll > 12)


    # isolation cuts
    outdict['l1_passes_tight_iso'] = data.relIso1 < 0.15
    outdict['l2_passes_tight_iso'] = data.relIso2 < 0.15
    # outdict['l1_passes_tight_iso'] = (data.relIso1 < 0.2)
    # outdict['l2_passes_tight_iso'] = (data.relIso2 < 0.2)
    outdict['l1_passes_loose_iso'] = (data.relIso1 > 0.2) & (data.relIso1 < 0.3)
    outdict['l2_passes_loose_iso'] = (data.relIso2 > 0.2) & (data.relIso2 < 0.3)
    outdict['l1_passes_ctrl_iso'] = outdict['l1_passes_loose_iso'] & ~outdict['l1_passes_tight_iso']
    outdict['l2_passes_ctrl_iso'] = outdict['l2_passes_loose_iso'] & ~outdict['l2_passes_tight_iso']

    outdict['isolation_sig'] = outdict['l1_passes_tight_iso'] & outdict['l2_passes_tight_iso'] & outdict['preselection']
    outdict['isolation_3l'] = outdict['l1_passes_tight_iso'] & outdict['l2_passes_tight_iso'] & outdict['3leptons']
    outdict['isolation_ss'] = outdict['l1_passes_tight_iso'] & outdict['l2_passes_tight_iso'] & outdict['same_sign']

    outdict['isolation_sig_ee'] = outdict['l1_passes_tight_iso'] & outdict['l2_passes_tight_iso'] & outdict['preselection'] & outdict['ee']
    outdict['isolation_sig_mumu'] = outdict['l1_passes_tight_iso'] & outdict['l2_passes_tight_iso'] & outdict['preselection'] & outdict['mumu']
    outdict['isolation_sig_emu'] = outdict['l1_passes_tight_iso'] & outdict['l2_passes_tight_iso'] & outdict['preselection'] & outdict['emu']

    outdict['isolation_ctrl'] = ( (outdict['l1_passes_tight_iso'] & outdict['l2_passes_ctrl_iso']) | (outdict['l1_passes_ctrl_iso'] & outdict['l2_passes_tight_iso']) )\
                                  & outdict['same_sign']

    outdict['sf'] = outdict['ee'] | outdict['mumu']
    outdict['of'] = outdict['emu']


    outdict['z_window'] = (data.mll < 106) & (data.mll > 76) & (abs(data.pdg1) == abs(data.pdg2))
    outdict['off_z_window'] = (((data.mll > 106) | (data.mll < 76)) | (abs(data.pdg1) != abs(data.pdg2))) & (data.mll > 12)
    # outdict['off_z_window'] = ~outdict['z_window'] & (data.mll > 12)

    outdict['pass_met'] = data.metPt > 60
    outdict['fail_met'] = data.metPt < 60

    # b-tag cuts
    outdict['0_bjets'] = data.nbjets == 0
    outdict['1_bjets'] = data.nbjets == 1
    outdict['2_bjets'] = data.nbjets == 2
    outdict['more_bjets'] = data.nbjets > 2
    outdict['bjets_sig'] = data.nbjets == 0
    outdict['bjets_ctrl'] = data.nbjets >= 1

    # combine to get signal and control samples
    # note we still haven't applied the mct cut
    outdict['wjets_ctrl'] = outdict['bjets_sig'] & outdict['isolation_ctrl'] & outdict['pass_met'] & outdict['off_z_window']
    outdict['wjets_ss_ctrl'] = outdict['bjets_sig'] & outdict['isolation_ss'] & outdict['pass_met'] & outdict['off_z_window']

    outdict['1tag_ctrl'] = outdict['1_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window']
    outdict['2tag_ctrl'] = outdict['2_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window']
    outdict['top_ctrl'] = outdict['bjets_ctrl'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window']
    outdict['moretag_ctrl'] = outdict['more_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window']
    outdict['z_ctrl'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['z_window']
    outdict['wz_ctrl'] = outdict['bjets_sig'] & outdict['isolation_3l'] & outdict['pass_met']
    outdict['sig'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window']

    outdict['wjets_ctrl_sf'] = outdict['bjets_sig'] & outdict['isolation_ctrl'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['sf']
    outdict['1tag_ctrl_sf'] = outdict['1_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['sf']
    outdict['2tag_ctrl_sf'] = outdict['2_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['sf']
    outdict['top_ctrl_sf'] = outdict['bjets_ctrl'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['sf']
    outdict['moretag_ctrl_sf'] = outdict['more_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['sf']
    outdict['z_ctrl_sf'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['z_window'] & outdict['sf']
    outdict['z_ctrl_0met'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['z_window'] & outdict['sf']
    outdict['z_ctrl_30met'] = outdict['bjets_sig'] & outdict['isolation_sig'] & (data.metPt > 30) & outdict['z_window'] & outdict['sf']
    outdict['sig_sf'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['sf']

    outdict['wjets_ctrl_of'] = outdict['bjets_sig'] & outdict['isolation_ctrl'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['of']
    outdict['1tag_ctrl_of'] = outdict['1_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['of']
    outdict['2tag_ctrl_of'] = outdict['2_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['of']
    outdict['top_ctrl_of'] = outdict['bjets_ctrl'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['of']
    outdict['moretag_ctrl_of'] = outdict['more_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['of']
    outdict['z_ctrl_of'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['z_window'] & outdict['of']
    outdict['sig_of'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['of']

    outdict['wjets_ctrl_ee'] = outdict['bjets_sig'] & outdict['isolation_ctrl'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['ee']
    outdict['1tag_ctrl_ee'] = outdict['1_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['ee']
    outdict['2tag_ctrl_ee'] = outdict['2_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['ee']
    outdict['top_ctrl_ee'] = outdict['bjets_ctrl'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['ee']
    outdict['z_ctrl_ee'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['z_window'] & outdict['ee']
    outdict['sig_ee'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['ee']

    outdict['wjets_ctrl_mumu'] = outdict['bjets_sig'] & outdict['isolation_ctrl'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['mumu']
    outdict['1tag_ctrl_mumu'] = outdict['1_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['mumu']
    outdict['2tag_ctrl_mumu'] = outdict['2_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['mumu']
    outdict['top_ctrl_mumu'] = outdict['bjets_ctrl'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['mumu']
    outdict['z_ctrl_mumu'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['z_window'] & outdict['mumu']
    outdict['sig_mumu'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['mumu']

    outdict['wjets_ctrl_emu'] = outdict['bjets_sig'] & outdict['isolation_ctrl'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['emu']
    outdict['1tag_ctrl_emu'] = outdict['1_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['emu']
    outdict['2tag_ctrl_emu'] = outdict['2_bjets'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['emu']
    outdict['top_ctrl_emu'] = outdict['bjets_ctrl'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['emu']
    outdict['z_ctrl_emu'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['z_window'] & outdict['emu']
    outdict['sig_emu'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met'] & outdict['off_z_window'] & outdict['emu']

    outdict['mct_low'] = (data.mctperp > 10.)
    outdict['mct_high'] = data.mctperp > mctcut

    outdict['sig_mct_low'] = outdict['sig'] & outdict['mct_low']
    outdict['sig_mct_high'] = outdict['sig'] & outdict['mct_high']
    outdict['wjets_mct_low'] = outdict['wjets_ctrl'] & outdict['mct_low']
    outdict['wjets_mct_high'] = outdict['wjets_ctrl'] & outdict['mct_high']
    outdict['top_mct_low'] = outdict['top_ctrl'] & outdict['mct_low']
    outdict['top_mct_high'] = outdict['top_ctrl'] & outdict['mct_high']
    outdict['1tag_mct_low'] = outdict['1tag_ctrl'] & outdict['mct_low']
    outdict['2tag_mct_low'] = outdict ['2tag_ctrl'] & outdict['mct_low']
    outdict['z_mct_low'] = outdict['z_ctrl'] & outdict['mct_low']
    outdict['z_mct_high'] = outdict['z_ctrl'] & outdict['mct_high']

    outdict['sig_mct_low_sf'] = outdict['sig'] & outdict['mct_low'] & outdict['sf']
    outdict['sig_mct_high_sf'] = outdict['sig'] & outdict['mct_high'] & outdict['sf']
    outdict['wjets_mct_low_sf'] = outdict['wjets_ctrl'] & outdict['mct_low'] & outdict['sf']
    outdict['wjets_mct_high_sf'] = outdict['wjets_ctrl'] & outdict['mct_high'] & outdict['sf']
    outdict['top_mct_low_sf'] = outdict['top_ctrl'] & outdict['mct_low'] & outdict['sf']
    outdict['top_mct_high_sf'] = outdict['top_ctrl'] & outdict['mct_high'] & outdict['sf']
    outdict['1tag_mct_low_sf'] = outdict['1tag_ctrl'] & outdict['mct_low'] & outdict['sf']
    outdict['2tag_mct_low_sf'] = outdict ['2tag_ctrl'] & outdict['mct_low'] & outdict['sf']
    outdict['z_mct_low_sf'] = outdict['z_ctrl'] & outdict['mct_low'] & outdict['sf']
    outdict['z_mct_high_sf'] = outdict['z_ctrl'] & outdict['mct_high'] & outdict['sf']

    outdict['sig_mct_low_of'] = outdict['sig'] & outdict['mct_low'] & outdict['of']
    outdict['sig_mct_high_of'] = outdict['sig'] & outdict['mct_high'] & outdict['of']
    outdict['wjets_mct_low_of'] = outdict['wjets_ctrl'] & outdict['mct_low'] & outdict['of']
    outdict['wjets_mct_high_of'] = outdict['wjets_ctrl'] & outdict['mct_high'] & outdict['of']
    outdict['top_mct_low_of'] = outdict['top_ctrl'] & outdict['mct_low'] & outdict['of']
    outdict['top_mct_high_of'] = outdict['top_ctrl'] & outdict['mct_high'] & outdict['of']
    outdict['1tag_mct_low_of'] = outdict['1tag_ctrl'] & outdict['mct_low'] & outdict['of']
    outdict['2tag_mct_low_of'] = outdict ['2tag_ctrl'] & outdict['mct_low'] & outdict['of']
    outdict['z_mct_low_of'] = outdict['z_ctrl'] & outdict['mct_low'] & outdict['of']
    outdict['z_mct_high_of'] = outdict['z_ctrl'] & outdict['mct_high'] & outdict['of']

    try:
        outdict['pass_met_scaleup'] = data.metPt_up > 60
        outdict['pass_met_scaledown'] = data.metPt_down > 60

        outdict['sig_scaleup'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met_scaleup'] & outdict['off_z_window']
        outdict['sig_scaledown'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met_scaledown'] & outdict['off_z_window']

        outdict['sig_scaleup_sf'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met_scaleup'] & outdict['off_z_window'] & outdict['sf']
        outdict['sig_scaledown_sf'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met_scaledown'] & outdict['off_z_window'] & outdict['sf']

        outdict['sig_scaleup_of'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met_scaleup'] & outdict['off_z_window'] & outdict['of']
        outdict['sig_scaledown_of'] = outdict['bjets_sig'] & outdict['isolation_sig'] & outdict['pass_met_scaledown'] & outdict['off_z_window'] & outdict['of']

        outdict['mct_high_scaleup'] = data.mctperp_up > mctcut
        outdict['mct_high_scaledown'] = data.mctperp_down > mctcut
        outdict['sig_mct_high_scaleup'] = outdict['sig_scaleup'] & outdict['mct_high_scaleup']
        outdict['sig_mct_high_scaledown'] = outdict['sig_scaledown'] & outdict['mct_high_scaledown']
        outdict['sig_mct_high_scaleup_sf'] = outdict['sig_scaleup'] & outdict['mct_high_scaleup'] & outdict['sf']
        outdict['sig_mct_high_scaledown_sf'] = outdict['sig_scaledown'] & outdict['mct_high_scaledown'] & outdict['sf']
        outdict['sig_mct_high_scaleup_of'] = outdict['sig_scaleup'] & outdict['mct_high_scaleup'] & outdict['of']
        outdict['sig_mct_high_scaledown_of'] = outdict['sig_scaledown'] & outdict['mct_high_scaledown'] & outdict['of']
    except AttributeError:
        pass

    return outdict
