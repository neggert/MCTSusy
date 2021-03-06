from selection import *

def get_mc_truth( data, mctcut, flavor) :
    sel = get_samples( data, mctcut)
    truth_high_dict = {}
    truth_high_dict['Total'] = sum(data[sel['sig_mct_high'+flavor]].weight)
    truth_high_dict['Top'] = sum(data[sel['sig_mct_high'+flavor] & (data.mc_cat == 'top')].weight)
    truth_high_dict['VV'] = sum(data[sel['sig_mct_high'+flavor] & ((data.mc_cat=="WV") | (data.mc_cat=="ZZ"))].weight)
    truth_high_dict['WW'] = sum(data[sel['sig_mct_high'+flavor] & (data.mctype.isin(['WWW','WWG','WWTo2L2Nu']))].weight)
    truth_high_dict['WZ'] = sum(data[sel['sig_mct_high'+flavor] & (data.mctype=='WZTo3LNu')].weight)
    truth_high_dict['ZZ'] = sum(data[sel['sig_mct_high'+flavor] & (data.mctype=='ZZTo2L2Nu')].weight)
    truth_high_dict['DY'] = sum(data[sel['sig_mct_high'+flavor] & (data.mc_cat=="DY")].weight)
    truth_high_dict['W'] = sum(data[sel['sig_mct_high'+flavor] & (data.mc_cat=="fake")].weight)


    truth_low_dict = {}
    truth_low_dict['Total'] = sum(data[sel['sig_mct_low'+flavor]].weight)
    truth_low_dict['Top'] = sum(data[sel['sig_mct_low'+flavor] & (data.mc_cat == 'top')].weight)
    truth_low_dict['VV'] = sum(data[sel['sig_mct_low'+flavor] & ((data.mc_cat=="WV") | (data.mc_cat=="ZZ"))].weight)
    truth_low_dict['WW'] = sum(data[sel['sig_mct_low'+flavor] & (data.mctype.isin(['WWW','WWG','WWTo2L2Nu']))].weight)
    truth_low_dict['WZ'] = sum(data[sel['sig_mct_low'+flavor] & (data.mctype=='WZTo3LNu')].weight)
    truth_low_dict['ZZ'] = sum(data[sel['sig_mct_low'+flavor] & (data.mctype=='ZZTo2L2Nu')].weight)
    truth_low_dict['DY'] = sum(data[sel['sig_mct_low'+flavor] & (data.mc_cat=="DY")].weight)
    truth_low_dict['W'] = sum(data[sel['sig_mct_low'+flavor] & (data.mc_cat=="fake")].weight)

    return {'low':truth_low_dict, 'high':truth_high_dict}
