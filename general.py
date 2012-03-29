import pylab

from data_handling import *

def get_dist_mct( hdf_file, isolated=True ) :
    if isolated :
        iterator = load_data_iterator_cut(hdf_file, """relIso2 < 0.17""")
    else :
        iterator = load_data_iterator( hdf_file )
    results = [data["mct"] for data in iterator]
    return results

def draw_mct_hist( hd5_file ) :
    """Draw mct and relIso hists from files"""
    data_mct = get_dist_mct( hd5_file )
    pylab.figure()
    pylab.hist(data_mct, bins=20)
    pylab.xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")

def get_dist_relIso( hdf_file ) :
    """Return a list of the relIso of the second most isolated lepton for each event in files"""
    return [data["relIso2"] for data in load_data_iterator(hdf_file)]

def draw_relIso_hist( hd5_file ) :
    data_relIso = get_dist_relIso(hd5_file)
    pylab.figure()
    pylab.hist(data_relIso, bins=20, range=(0, 2))
    pylab.xlabel("relIso")
    pylab.show()

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

def make_numpy_hist( data, weight) :
    return pylab.histogram(data, bins=20, range=(1,280), weights=[weight for d in data])


