from pandas import HDFStore
from selection import get_samples
import sys

sf_only = False

if len(sys.argv) >= 4:
    sf_only = bool(sys.argv[3])

s = HDFStore(sys.argv[1])
d = s['data']

sd = get_samples(d)

dsig_of = d[sd['sig_mct_low_of']]

masses_of = set()
for g, _ in dsig_of.groupby(['mass1', 'mass2']):
	masses_of.add(tuple(map(int, g)))

dsig_sf = d[sd['sig_mct_low_sf']]

masses_sf = set()
for g, _ in dsig_sf.groupby(['mass1', 'mass2']):
    masses_sf.add(tuple(map(int, g)))

if sf_only:
    masses = masses_sf
else:
    masses = masses_sf.intersection(masses_of)

if sf_only:
    masses = filter( lambda x: x[0] - x[1] > 50, masses)
else:
    masses = filter( lambda x: x[0] - x[1] > 75, masses)

if sf_only:
    masses = filter( lambda x: x[0] <= 475, masses)

print " ".join(map(lambda x: "limits/"+sys.argv[2]+"_"+str(int(x[0]))+"_"+str(int(x[1]))+"_combined_meas_model.root", masses))

s.close()
