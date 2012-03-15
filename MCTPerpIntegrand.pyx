from libc.math cimport sin, exp, log, sqrt, fabs

import numpy
cdef double pi = numpy.pi

cdef class MCTPerp_integrand(object):
    """Integrand for calculating mct perp likelihood"""
    cdef double ptprod
    cdef double delta
    cdef double phiUp
    cdef double sPhiUp
    cdef double frac
    cdef double ep
    def __init__( self ):
        """Empty constructor"""

    def set_data (self, data, double frac, double endpoint) :
        self.ptprod = data["pt1"]*data["pt2"]
        self.delta = data["delta"]
        self.phiUp = data["phiUp"]
        self.sPhiUp = data["sPhiUp"]
        self.frac = frac
        self.ep = endpoint
    
    cdef double pdf( self, double mct) :
        """return mct pdf for mct(phi)"""
        if (mct > 0. ) : return (1-self.frac)*4*mct/self.ep**2*log(self.ep/mct)
        else : return 0.
    
    cdef double response( self, double phi ) :
        """normal distribution centered at phiUp evaluated at phi"""
        cdef double dPhi = fabs((phi-self.phiUp))
        if dPhi > pi : dPhi = 2*pi-dPhi
        return 1./sqrt(2*pi)/self.sPhiUp*exp(-dPhi*dPhi/2/self.sPhiUp/self.sPhiUp)
    
    cdef double jacobian( self, double phi, double mct) :
        """Return the jacobian factor dMCT/dPhi evaluated at phi"""
        return 2*fabs(self.ptprod*sin(2*phi))/mct
    
    cdef double mct( self, double phi ) :
        """calculate mct for this event at the given phi. Decorated by a lru cache because
        this calculation will be performed several times for the same input value"""
        return 2*sqrt(self.ptprod*max([sin(phi-self.delta)*sin(phi+self.delta), 0.]))

    def __call__(self, phi) :
        """Putting it all together"""
        cdef double mct = self.mct( phi )
        if mct > self.ep : return 0.
        elif mct > 0. : return self.pdf(mct)*self.response(phi)*self.jacobian(phi, mct)
        else : return self.frac/4/self.delta*self.response(phi)