import prep_hists
import bkg_fit

from pandas import *

def find_bootstrap_dist(filename="bootstrap.root", num_pes=1000):
    outputs = []
    for i in xrange(num_pes):
        prep_hists.create_data_file(filename, 14, (10,150), i+1)

        r = bkg_fit.run_bonly_fit("sig_chi.root", 200, 0, ['of', 'sf'], ncpu=8, get_p=False, data_prefix="data_bs_"+str(i+1), data_file_name=filename)
        outputs.append(r)

    return outputs

def process_results(res):
    of_data = [zip(*r['of'].values())[0] for r in res]
    sf_data = [zip(*r['sf'].values())[0] for r in res]

    df_of = DataFrame(of_data)
    df_sf = DataFrame(sf_data)

    df_of.columns = Index([name+'_of' for name in r['of'].keys()])
    df_sf.columns = Index([name+'_sf' for name in r['sf'].keys()])

    df = df_of.join(df_sf)

    return df

