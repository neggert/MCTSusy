import pylab

"""
Functions useful for cut and count analysis
"""

def get_events_pass_mct (mct_data, weight, mct_cut) :
    """Return the number of events with mct > mct_cut"""
    n_events = len([d for d in mct_data if d > mct_cut])
    return n_events*weight

def get_events_vs_mct_cut (hd5_file, weight) :
    """Return a tuple of the number of events passing vs mct cut"""
    data_mct = get_dist_mct( hd5_file )

    cuts = range(60, 120)
    passed = []
    for cut in cuts :
        passed.append(get_events_pass_mct(data_mct, weight, cut))
    return (cuts, passed)

def draw_num_events_past_mct_cut( hd5_file, weight ) :
    """Draw number of events passing cut vs cut value for files"""
    data = get_events_vs_mct_cut( hd5_file, weight )
    pylab.scatter(*data)
    pylab.ylim = (0, 20)
    pylab.xlabel("MCT cut")
    pylab.ylabel("Number of events")
    return data

def get_sig_over_sqrt_bkg_vs_mct_cut ( signal_tpl, bkg_tpl_list ) :
    """Get significance = signal/sqrt(bkg) vs mct cut, bkg is an list of tuples,
    each tuple contains a background pickle file and it's weight"""
    cut, s = pylab.array(get_events_vs_mct_cut( *signal_tpl ))
    b = pylab.zeros(len(s))
    for bkg_tpl in bkg_tpl_list :
        b += pylab.array(get_events_vs_mct_cut( *bkg_tpl ))[1]
    print b
    return cut, s/pylab.sqrt(b)

def passes_iso_cut ( event ) :
    for i in [1,2] :
        if abs(event["pdg"+str(i)]) == 11 :
            if event['relIso'+str(i)] > 0.17 :
                return False
        elif abs(event["pdg"+str(i)] == 13) :
            if event['relIso'+str(i)] > 0.20 :
                return False
    return True

def passes_inverted_iso_cut ( event ) :
    """ One lepton passes isolation cut, the other fails"""
    npass = 0
    nfail = 0
    for i in [1,2] :
        if abs(event["pdg"+str(i)]) == 11 :
            if 0.17 < event['relIso'+str(i)] < 0.4 :
                nfail += 1
            elif event['relIso'+str(i)] < 0.17 :
                npass += 1
        elif abs(event["pdg"+str(i)] == 13) :
            if 0.20 < event['relIso'+str(i)] < 0.4 :
                nfail += 1
            elif event['relIso'+str(i)] < 0.20 :
                npass += 1
    return (npass==1 and nfail==1)

def passes_mct_cut( event, cut=100.) :
    return event['mct'] > cut

def passes_bjet_veto( event ) :
    return event['nbjets'] == 0

def passes_bjet_cut( event, ntags = 2 ) :
    return event['nbjets'] >= ntags

def passes_sig_cuts( event, cut=100. ) :
    return passes_mct_cut( event, cut) and passes_iso_cut( event ) and passes_bjet_veto( event )

def passes_wjets_control_cuts( event, cut=100. ) :
    return passes_mct_cut( event, cut) and passes_inverted_iso_cut( event ) and passes_bjet_veto( event )

def passes_ttbar_control_cuts( event, cut=100., ntags=2) :
    return passes_mct_cut( event, cut) and passes_iso_cut( event ) and passes_bjet_cut( event, ntags)

def count_events( data_it ) :
    """Count events in the iterator. Warning! Destroys the iterator!"""
    return sum(1 for _ in data_it)