from config.data import *
import numpy as np

def im(pta, phia, etaa, ptb, phib, etab):
	return np.sqrt(2*pta*ptb*(np.cosh(etaa-etab)-np.cos(phia-phib)))


mc['m23'] = im(mc.pt2, mc.phi2, mc.eta2, mc.pt3, mc.phi3, mc.eta3)
mc['m13'] = im(mc.pt1, mc.phi1, mc.eta1, mc.pt3, mc.phi3, mc.eta3)
data['m23'] = im(data.pt2, data.phi2, data.eta2, data.pt3, data.phi3, data.eta3)
data['m13'] = im(data.pt1, data.phi1, data.eta1, data.pt3, data.phi3, data.eta3)