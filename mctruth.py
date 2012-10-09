from selection import *

def get_mc_truth( data, mctcut) :
    sel = get_samples( data, mctcut)
    truth_high_dict = {}
    truth_high_dict['Total'] = sum(data[sel['sig_mct_high']].weight)
    truth_high_dict['Top'] = sum(data[sel['sig_mct_high'] & (data.mc_cat == 'top')].weight)
    truth_high_dict['WV'] = sum(data[sel['sig_mct_high'] & ((data.mc_cat=="WV") | (data.mc_cat=="ZZ"))].weight)
    truth_high_dict['DY'] = sum(data[sel['sig_mct_high'] & (data.mc_cat=="DY")].weight)
    truth_high_dict['W'] = sum(data[sel['sig_mct_high'] & (data.mc_cat=="fake")].weight)


    truth_low_dict = {}
    truth_low_dict['Total'] = sum(data[sel['sig_mct_low']].weight)
    truth_low_dict['Top'] = sum(data[sel['sig_mct_low'] & (data.mc_cat == 'top')].weight)
    truth_low_dict['WV'] = sum(data[sel['sig_mct_low'] & ((data.mc_cat=="WV") | (data.mc_cat=="ZZ"))].weight)
    truth_low_dict['DY'] = sum(data[sel['sig_mct_low'] & (data.mc_cat=="DY")].weight)
    truth_low_dict['W'] = sum(data[sel['sig_mct_low'] & (data.mc_cat=="fake")].weight)

    return {'low':truth_low_dict, 'high':truth_high_dict}
