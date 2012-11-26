import matplotlib
# matplotlib.use("pdf")
# matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['text.latex.preamble'] = '\usepackage{mathpazo}, \usepackage{upgreek}'
# matplotlib.rcParams['font.size'] = 14
# matplotlib.rcParams['figure.figsize'] = 8, 8
from matplotlib.pylab import *


from pandas import *
import sys
sys.path.append("../")
from CMSPyLibs.plot import *
from selection import *
from stackplot import check_for_data

store = HDFStore("Data/mc2012.hdf5")

smc = get_samples(store['all_cat'])

selection = smc['sig']

variable = "mctperp"
bins = 9
plotrange = (40, 400)

bkgtpl = []
bkgwtpl = []
bkgltpl = []

groups = store['all_cat'][selection].groupby('mc_cat')

group_order = ['top', 'WV', 'ZZ', 'DY', 'fake']

for name in group_order
    bkgtpl.append( store['all_cat'][selection][groups.groups()[name]][variable])
    bkgwtpl.append( store['all_cat'][selection][groups.groups()[name]].weight)
    bkgltpl.append(name)

print "Plotting", len(bkgtpl), "backgrounds..."
figure()
fig = subplot(111)
fig.set_yscale('log', nonposy='clip')
hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl)

xlabel("$M_{\mathrm{CT}\perp}$ (GeV)")
legend()
show()
savefig("huh.pdf")

store.close()
