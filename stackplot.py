from matplotlib.pyplot import hist

def stackplot( data, selection, variable, bins, range) :
    """Make a stackhist in the given variable"""
    bkgtpl = (data[selection & (data.mctype=="WW")][variable], data[selection & (data.mctype=="ttbar")][variable],\
        data[selection & (data.mctype=="tW")][variable], data[selection & (data.mctype=="DY")][variable],\
        data[selection & (data.mctype=="wjets")][variable])
    bkgwtpl = (data[selection & (data.mctype=="WW")].weight, data[selection & (data.mctype=="ttbar")].weight,\
        data[selection & (data.mctype=="tW")].weight, data[selection & (data.mctype=="DY")].weight,\
        data[selection & (data.mctype=="wjets")].weight)
    bkgltpl = ("WW", "t\=t", "tW", "Z/$\gamma^*$", "W+Jets" )
    hist(bkgtpl, weights=bkgwtpl, histtype="stepfilled", stacked=True, rwidth=1, bins=bins, range=range, label=bkgltpl)
    hist( data[selection & (data.mctype=="ChiX0W_300_100")][variable], weights=data[selection & (data.mctype=="ChiX0W_300_100")].weight, \
        histtype="step", fill=False, bins=bins, range=range, label=r"$\chi^+\chi^- \to \chi^0\mathrm{W}^+\chi^0\mathrm{W}^-$", lw=2, ec="k")
    hist( data[selection & (data.mctype=="ChilSnu_300_100")][variable], weights=data[selection & (data.mctype=="ChilSnu_300_100")].weight, \
        histtype="step", fill=False, bins=bins, range=range, label=r"$\chi^+\chi^- \to \tilde{\nu}\ell^+\tilde{\nu}\ell^-$", lw=2, ec="#aaaaaa")