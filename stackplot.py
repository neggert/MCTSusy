from matplotlib.pyplot import hist
from CMSPyLibs.plot import *
from copy import copy


def check_for_data(data, selection, mctype, variable, plotrange) :
    out =  data[selection & (data.mctype.isin(mctype)) & (data[variable] >= plotrange[0]) & (data[variable] <= plotrange[1])][variable].count() > 0
    if not out :
        print "Not plotting", mctype
    return out

def stackplot( data, real_data, selection, real_selection, variable, bins, plotrange) :
    """Make a stackhist in the given variable"""
    dysamples = ['DYToEE_M10', 'DYToEE_M20', 'DYToMuMu_M10', 'DYToMuMu_M20', 'DYToTauTau_M10', 'DYToTauTau_M20', 'DYJets_M50']
    bkgtpl = []
    bkgwtpl = []
    bkgltpl = []
    bkg = "WW"
    if check_for_data( data, selection, [bkg,], variable, plotrange) :
        bkgtpl.append( data[selection & (data.mctype==bkg)][variable])
        bkgwtpl.append( data[selection & (data.mctype==bkg)].weight)
        bkgltpl.append("WW")
    bkg = "WZ"
    if check_for_data( data, selection, [bkg,], variable, plotrange) :
        bkgtpl.append( data[selection & (data.mctype==bkg)][variable])
        bkgwtpl.append( data[selection & (data.mctype==bkg)].weight)
        bkgltpl.append("WZ")
    bkg = "ZZ"
    if check_for_data( data, selection, [bkg,], variable, plotrange) :
        bkgtpl.append( data[selection & (data.mctype==bkg)][variable])
        bkgwtpl.append( data[selection & (data.mctype==bkg)].weight)
        bkgltpl.append("ZZ")
    bkg = "ttbar"
    if check_for_data( data, selection, [bkg,], variable, plotrange) :
        bkgtpl.append( data[selection & (data.mctype==bkg)][variable])
        bkgwtpl.append( data[selection & (data.mctype==bkg)].weight)
        bkgltpl.append("t\=t")
    bkg = "tW"
    if check_for_data( data, selection, [bkg,], variable, plotrange) :
        bkgtpl.append( data[selection & (data.mctype==bkg)][variable])
        bkgwtpl.append( data[selection & (data.mctype==bkg)].weight)
        bkgltpl.append("tW")
    bkg = "DY"
    if check_for_data( data, selection, dysamples, variable, plotrange) :
        bkgtpl.append( data[selection & data.mctype.isin(dysamples)][variable])
        bkgwtpl.append( data[selection & data.mctype.isin(dysamples)].weight)
        bkgltpl.append("Z/$\gamma^*$")
    bkg = "wjets"
    if check_for_data( data, selection, [bkg,], variable, plotrange) :
        bkgtpl.append( data[selection & (data.mctype==bkg)][variable])
        bkgwtpl.append( data[selection & (data.mctype==bkg)].weight)
        bkgltpl.append("W+Jets")
    print "Plotting", len(bkgtpl), "backgrounds..."
    hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=plotrange, label=bkgltpl)
    #hist( data[selection & (data.mctype=="300_100")][variable], weights=data[selection & (data.mctype=="300_100")].weight, \
    #    histtype="step", fill=False, bins=bins, plotrange=plotrange, label=r"$\chi^+\chi^- \to \chi^0\mathrm{W}^+\chi^0\mathrm{W}^-$", lw=2, ec="k")
    he = hist_errorbars( real_data[real_selection][variable], xerrs=False, bins=bins, range=plotrange)
    he.set_label("Data")

