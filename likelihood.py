
import pylab
from pylab import pi, sqrt, sin, log, exp
from scipy.integrate import quad
from general import *
import event_counter
import events.dilepton_event
import tables
import itertools
import pdb
from MCTPerpIntegrand import MCTPerp_integrand

class LikelihoodCalculator(object):
    """Calculator for the likelihood of the signal fraction given the data"""
    sm_endpoint = 80.4
    integrand = MCTPerp_integrand()
    def __init__(self, hd5_files, weights):
        """Intialize with a list of hd5_files"""
        self.files = [tables.openFile(hd5_file, "r") for hd5_file in hd5_files]
        self.table = [f.root.fit.fitdata for f in self.files]
        self.rand_state = pylab.get_state() # save this so we can produce the same random numbers every time
        self.weights = weights

        self.nevents = len([x["mctype"] for x in self.get_data_iterator()])
    
    def sampled_iter(self, sequence, p) :
        for elem in sequence:
            if pylab.rand() < p :
                yield elem


    def get_data_iterator(self ) :
        """Make a copy of the iterator"""
        iterators = [t.where("""relIso2 < 0.17""") for t in self.table]
        pylab.set_state(self.rand_state) # Make sure each time this is run it produces the same random numbers
        filtered_iterators = [self.sampled_iter(i, w) for w,i in zip(self.weights, iterators)]
        return itertools.chain(*filtered_iterators)

    
    def ttbar_likelihood( self, data) :
        """return the likelihood for a single event from a single distribution"""
        delta = data["delta"]
        self.integrand.set_data( data, 0.5, self.sm_endpoint)
        points = [-pi+delta, -delta, delta, pi-delta]
        result = quad( self.integrand, -pi, pi, limit=100, points=points)
        if result == 0. :
            pdb.set_trace()
        return result[0]

    def get_nll( self ) :
        """
        return the likelihood from two summed distributions with endpoint 1 and 2.
        The first distribution has fraction (1-s) and the second has fraction s.
        Note that this implementation is probably slow, since we have to call
        single_dist_likelihood twice, and thus loop through all events twice.
        """
        ll =0.
        for data in self.get_data_iterator() :
            l = log(self.ttbar_likelihood( data ))
            ll -= l
        return ll
    
    def __eval__ (self, s) :
        """
        Return the marginalized likelihood of a signal value s for endpoint1 = sm_endpoint
        """
        pass
    

def run_fit( data_tuples ) :
    """tuple should should be of the form (pickle_file, weight)"""
    pass
