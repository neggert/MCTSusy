import mt2.mt2_bisect as mt2
import numpy as np


def mt2calc(event):
    m = mt2.mt2()

    m.set_momenta((0, event['pt1']*np.sin(event['phi1']), event['pt1']*np.cos(event['phi1'])),
                  (0, event['pt2']*np.sin(event['phi2']), event['pt2']*np.cos(event['phi2'])),
                  (0, event['metPt']*np.sin(event['metPhi']), event['metPt']*np.cos(event['metPhi'])))
    m.set_mn(0.)

    return m.get_mt2()
