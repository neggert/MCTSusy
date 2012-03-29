import pylab

from data_handling import *

def get_dist_mct( hdf_file, isolated=True ) :
    if isolated :
        iterator = load_data_iterator_iso(hdf_file)
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

def make_numpy_hist( data, weight) :
    return pylab.histogram(data, bins=20, range=(1,280), weights=[weight for d in data])


