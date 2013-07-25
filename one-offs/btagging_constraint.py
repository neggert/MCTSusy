"""
This script is meant to test my idea for finding the number of un-tagged
top events by comparing the results of two different b-taggers.
"""
from pandas import *
from math import log
import selection
reload(selection)

# load the dataset I made for this occasion
s = HDFStore("work/mc/b-tagging.hdf5")

mc = s['data']

smc = selection.get_samples(mc)

# get the data that passes all cuts but the b-veto
sig_mc = mc[smc['sig_no_b_veto']&smc['mct_low']]

# now we get the numbers of events with different taggers

n_csvl = {}
for n_tags, data in sig_mc.groupby("n_jets_csvl"):
    n_csvl[n_tags] = sum(data.weight)

n_csvm = {}
for n_tags, data in sig_mc.groupby("n_jets_csvm"):
    n_csvm[n_tags] = sum(data.weight)

n_jpl = {}
for n_tags, data in sig_mc.groupby("n_jets_jpl"):
    n_jpl[n_tags] = sum(data.weight)

n_jpm = {}
for n_tags, data in sig_mc.groupby("n_jets_jpm"):
    n_jpm[n_tags] = sum(data.weight)

# Get efficiencies
n_gen_bjets = sum(sig_mc.n_gen_bjets * sig_mc.weight)

n_csvm_tags = sum(sig_mc.n_jets_csvm * sig_mc.weight)
n_jpm_tags = sum(sig_mc.n_jets_jpm * sig_mc.weight)

eff_csvm = 1. * n_csvm_tags / n_gen_bjets
eff_jpm = 1. * n_jpm_tags / n_gen_bjets

print "Efficiencies - CSVM: {}\t JPM: {}".format(eff_csvm, eff_jpm)


# don't forget to use the scale factors in data

# can we solve the system of equations using python?
# from import.categorize import categorize
# cats = sig_mc.mctype.apply(categorize)
# sig_mc.insert(len(sig_mc.columns), 'mc_cat', cats)

from minuit2 import Minuit2 as Minuit

eps_a = eff_csvm
N1_a = n_csvm[1]
N2_a = n_csvm[2]

eps_b = eff_jpm
N1_b = n_jpm[1]
N2_b = n_jpm[2]

p0 = lambda eps, A, ftt: ((1-A)**2+2*A*(1-A)*(1-eps)+A**2*(1-eps)**2)*ftt+(1-eps*A)*(1-ftt)
p1 = lambda eps, A, ftt: eps*A*(1-ftt)+(2*A*(1-A)*eps+2*A**2*eps*(1-eps))*ftt
p2 = lambda eps, A, ftt: ftt*eps**2*A**2


def logL(eps, N1, N2, Nt, A, ftt):
    return Nt*log(Nt)-(Nt-N1-N2)*log(Nt-N1-N2)+(Nt-N1-N2)*log(p0(eps, A, ftt))+N1*log(p1(eps, A, ftt))+N2*log(p2(eps, A, ftt))


def NLL(Nt, A, ftt):
    return -logL(eps_a, N1_a, N2_a, Nt, A, ftt) -logL(eps_b, N1_b, N2_b, Nt, A, ftt)

m = Minuit(NLL, Nt=5*N2_a, A=0.9, ftt=0.9)
m.limits['Nt']=(max([N1_a+N2_a, N1_b+N2_b]), 100000000)
m.limits['ftt'] = (0, 1)
m.limits['A'] = (0, 1)
m.migrad()

print m.values
print m.values['Nt']-N1_a-N2_a
print m.errors

print sum(sig_mc[(sig_mc.mctype.apply(str.lower).isin(["tw", 'ttw', 'ttbar', 'ttg', 'ttww', 'ttz'])) & (sig_mc.n_jets_csvm==0)].weight)
