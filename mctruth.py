from selection import *

def get_mc_truth( data, mctcut) :
    sel = get_samples( data, mctcut)
    truth_high_dict = {}
    truth_high_dict['Total'] = sum(data[sel['sig_mct_high']].weight)
    truth_high_dict['Top'] = sum(data[sel['sig_mct_high'] & ((data.mctype=="ttbar") | (data.mctype=='tW'))].weight)
    truth_high_dict['WV'] = sum(data[sel['sig_mct_high'] & ((data.mctype=="WW") | (data.mctype=="WZ"))].weight)
    truth_high_dict['ZZ'] = sum(data[sel['sig_mct_high'] & (data.mctype=="ZZ")].weight)
    truth_high_dict['DY'] = sum(data[sel['sig_mct_high'] & (data.mctype=="DY")].weight)
    truth_high_dict['W'] = sum(data[sel['sig_mct_high'] & (data.mctype=="DY")].weight)


    truth_low_dict = {}
    truth_low_dict['Total'] = sum(data[sel['sig_mct_low']].weight)
    truth_low_dict['Top'] = sum(data[sel['sig_mct_low'] & ((data.mctype=="ttbar") | (data.mctype=='tW'))].weight)
    truth_low_dict['WV'] = sum(data[sel['sig_mct_low'] & ((data.mctype=="WW") | (data.mctype=="WZ"))].weight)
    truth_low_dict['ZZ'] = sum(data[sel['sig_mct_low'] & (data.mctype=="ZZ")].weight)
    truth_low_dict['DY'] = sum(data[sel['sig_mct_low'] & (data.mctype=="DY")].weight)
    truth_low_dict['W'] = sum(data[sel['sig_mct_low'] & (data.mctype=="DY")].weight)

    return {'low':truth_low_dict, 'high':truth_high_dict}
