from pandas import *

from selection import *


s = HDFStore("work/mc/tchipm_gen_pt.hdf5")

d = s['data']

sd = get_samples(d)