import pylab
import numpy.random
from numpy import pi, sqrt, sin, log, exp, vectorize, array
from scipy.integrate import quad
import tables
import itertools
import pdb
from MCTPerpIntegrand import MCTPerp_integrand

class LikelihoodCalculator(object):
    """
    Calculator for the likelihood of the data. There is a random element to which
    data points are chosen, so results will vary from instance to instance, even
    with the same constructor. The seed is saved, so results within an instance will
    be consistant.
    """
    sm_endpoint = 80.4
    integrand = MCTPerp_integrand()
    def __init__(self, hd5_files, weights):
        """Intialize with a list of hd5_files"""
        self.files = [tables.openFile(hd5_file, "r") for hd5_file in hd5_files]
        self.tables = [f.root.fit.fitdata for f in self.files]
        self.weights = weights
        self.eventlists = self.resample()

        # self.nevents = len([x["mctype"] for x in self.get_data_iterator()])

    def resample(self) :
        """Sample event numbers with replacement"""
        eventlists = []
        for i, t in enumerate(self.tables) :
            eventlist = []
            # a somewhat roundabout way of getting the number of events that pass the relIso cut
            nevents_total = len([1 for d in t.read() if d['relIso2'] < 0.17])
            nevents_select = int(nevents_total*self.weights[i])
            for j in range(nevents_select) :
                eventnum = numpy.random.randint(nevents_total)
                eventlist.append(eventnum)
            eventlists.append(eventlist)
        return eventlists
    
    def sampled_iter(self, sequence, eventlist) :
        for event in eventlist:
            yield sequence[event]


    def get_data_iterator(self ) :
        """Make a copy of the iterator."""
        iterators = [[d for d in t.read() if d['relIso2'] < 0.17] for t in self.tables]
        filtered_iterators = [self.sampled_iter(i, l) for l,i in zip(self.eventlists, iterators)]
        return itertools.chain(*filtered_iterators)

    
    def ttbar_likelihood( self, data) :
        """return the likelihood for a single event to be from the ttbar distribution"""
        delta = data["delta"]
        self.integrand.set_data( data, 0.5, self.sm_endpoint)
        points = [-pi+delta, -delta, delta, pi-delta]
        result = quad( self.integrand, -pi, pi, limit=100, points=points)
        if result == 0. :
            pdb.set_trace()
        return result[0]

    def get_nll( self, cut=-1 ) :
        """
        Get the negative log likelihood summed over all events.
        """
        return sum(self.get_nlls(cut))

    def get_mcts( self ) :
        """
        return a list of mct values for each event in this sample
        """
        mcts = [data['mct'] for data in self.get_data_iterator()]
        return mcts

    def get_nlls( self, cut=-1 ) :
        """Return a list of nll values for each event in this sample, possibly with an MCT cut"""
        nlls = [-log(self.ttbar_likelihood( data )) for data in self.get_data_iterator() if data['mct'] > cut]
        return nlls
   
def get_210_likelihood( data, endpoint ) :
    """Get the MCTPerp 210 likelihood for an event given an endpoint"""
    delta = data["delta"]
    integrand = MCTPerp_integrand()
    integrand.set_data( data, 0.5, endpoint)
    points = [-pi+delta, -delta, delta, pi-delta]
    result = quad( integrand, -pi, pi, limit=100, points=points)[0]
    return result

def get_SplusB_NLL( data_list, s_endpoint, b_endpoint, N_s, N_b) :
    """(Extended) negative log likelihood for signal plus background"""
    s_likelihoods = array([get_210_likelihood(x, s_endpoint) for x in data_list])
    b_likelihoods = array([get_210_likelihood(x, b_endpoint) for x in data_list])
    return sum(-log(N_b*b_likelihoods+N_s*s_likelihoods)+N_s+N_b)

def get_B_NLL( data_list, b_endpoint, N_b) :
    """(Extended) negative log likelihood for signal plus background"""
    b_likelihoods = array([get_210_likelihood(x, b_endpoint) for x in data_list])
    return sum(-log(N_b*b_likelihoods)+N_b)
