import matplotlib.pyplot as plt
import numpy as np
import inspect

def hist_ratio( data_num, data_denom, weight_denom, xerrs=True, color="k", *args, **kwargs) :
    """Plot a histogram with error bars. Accepts any kwarg accepted by either numpy.histogram or pyplot.errorbar"""
    # pop off normed kwarg, since we want to handle it specially

    # retrieve the kwargs for numpy.histogram
    histkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(np.histogram).args :
            histkwargs[key] = value

    histvals1, binedges = np.histogram( data_num, **histkwargs )
    yerrs1 = np.sqrt(histvals1.tolist()) # no effing idea why tolist is necessary

    histvals2, binedges = np.histogram( data_denom, weights=weight_denom, **histkwargs )
    yerrs2 = np.sqrt(histvals2.tolist()) # no effing idea why tolist is necessary

    ratio = histvals1/histvals2
    ratio_err = np.sqrt((yerrs1/histvals2)**2+(histvals1/histvals2**2*yerrs2)**2)

    bincenters = (binedges[1:]+binedges[:-1])/2

    if xerrs :
        xerrs = (binedges[1]-binedges[0])/2
    else :
        xerrs = None

    # retrieve the kwargs for errorbar
    ebkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(plt.errorbar).args :
            histkwargs[key] = value
    out = plt.errorbar(bincenters, ratio, ratio_err, xerrs, fmt=".", color=color, **ebkwargs)

    if 'log' in kwargs.keys() :
        if kwargs['log'] :
            plt.yscale('log')

    if 'range' in kwargs.keys() :
        plt.xlim(*kwargs['range'])

    return out
