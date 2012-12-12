from IPython.parallel import Client

rc = Client()
lview = rc.load_balanced_view()
dview = rc[:]

with dview.sync_imports():
    import set_limits

import json

from config.data import chi, sel_chi



ncpu = 1
asymptotic = True
coarse = False
chans = ['sf', 'of']
sig_file = "sig_chi.root"

selected = chi[sel_chi['sig'] & (chi.mctperp>10.)]
groups = selected.groupby(['mass1', 'mass2'])

lims = {}

results = {}

for name, data in groups:
    m1, m2 = name
    m1 = int(m1)
    m2 = int(m2)

    results[(m1,m2)]= lview.apply_async(set_limits.run_limit, sig_file, m1, m2, chans, ncpu, asymptotic, coarse)

lview.wait(results.values())

for key in results.keys():
    res = results[key]
    lims[str(key)] = res.get()

f = open("chi_expected_limits.json", "w")

json.dump(lims, f, indent=2)

f.close()
