
def get_samples( data ) :
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

    # isolation cuts
    outdict['l1_passes_tight_iso'] = ((abs((data.pdg1.values) == 11) & (data.relIso1 < 0.17)) | ((abs(data.pdg1.values) == 13) & (data.relIso1 < 0.2)))
    outdict['l2_passes_tight_iso'] = ((abs((data.pdg2.values) == 11) & (data.relIso2 < 0.17)) | ((abs(data.pdg2.values) == 13) & (data.relIso2 < 0.2)))
    outdict['l1_passes_loose_iso'] = ((abs((data.pdg1.values) == 11) & (data.relIso1 < 0.4)) | ((abs(data.pdg1.values) == 13) & (data.relIso1 < 0.4)))
    outdict['l2_passes_loose_iso'] = ((abs((data.pdg2.values) == 11) & (data.relIso2 < 0.4)) | ((abs(data.pdg2.values) == 13) & (data.relIso2 < 0.4)))
    outdict['l1_passes_ctrl_iso'] = outdict['l1_passes_loose_iso'] & ~outdict['l1_passes_tight_iso']
    outdict['l2_passes_ctrl_iso'] = outdict['l2_passes_loose_iso'] & ~outdict['l2_passes_tight_iso']

    outdict['isolation_sig'] = outdict['l1_passes_tight_iso'] & outdict['l2_passes_tight_iso']
    outdict['isolation_ctrl'] = (outdict['l1_passes_tight_iso'] & outdict['l2_passes_ctrl_iso']) | (outdict['l1_passes_ctrl_iso'] & outdict['l2_passes_tight_iso'])

    # b-tag cuts
    outdict['0_bjets'] = data.nbjets == 0
    outdict['1_bjets'] = data.nbjets == 1
    outdict['2_bjets'] = data.nbjets == 2
    outdict['more_bjets'] = data.nbjets > 2
    outdict['bjets_sig'] = data.nbjets == 0
    outdict['bjets_ctrl'] = data.nbjets >= 1

    # combine to get signal and control samples
    # note we still haven't applied the mct cut
    outdict['wjets_ctrl'] = outdict['bjets_sig'] & outdict['isolation_ctrl']
    outdict['1tag_ctrl'] = outdict['1_bjets'] & outdict['isolation_sig']
    outdict['2tag_ctrl'] = outdict['2_bjets'] & outdict['isolation_sig']
    outdict['top_ctrl'] = outdict['bjets_ctrl'] & outdict['isolation_sig']
    outdict['iso_bjet_sig'] = outdict['bjets_sig'] & outdict['isolation_sig']

    outdict['mct_low'] = (data.mctperp > 5.) & (data.mctperp < 100.)
    outdict['mct_high'] = data.mctperp > 100.

    outdict['sig_mct_low'] = outdict['iso_bjet_sig'] & outdict['mct_low']
    outdict['sig_mct_high'] = outdict['iso_bjet_sig'] & outdict['mct_high']
    outdict['wjets_mct_low'] = outdict['wjets_ctrl'] & outdict['mct_low']
    outdict['wjets_mct_high'] = outdict['wjets_ctrl'] & outdict['mct_high']
    outdict['top_mct_low'] = outdict['top_ctrl'] & outdict['mct_low']
    outdict['top_mct_high'] = outdict['top_ctrl'] & outdict['mct_high']
    outdict['1tag_mct_low'] = outdict['1tag_ctrl'] & outdict['mct_low']
    outdict['2tag_mct_low'] = outdict ['2tag_ctrl'] & outdict['mct_low']

    return outdict