import sys
sys.path.append("..")
import BackgroundFit
import mctruth
reload(BackgroundFit)
reload(mctruth)
import numpy
import numpy.random
from collections import defaultdict
from data import *


def run_mc_closure(flavor):
    """
    Check MC closure for a range of MCT cuts. Writes a latex table containing the results
    """

    # all mc except diboson and DY
    test_mc = mc[((mc.mc_cat != 'WV') & (mc.mc_cat != 'ZZ') & (mc.mc_cat != 'DY'))]


    diboson_mc = mc[((mc.mc_cat == 'WV') | (mc.mc_cat == 'ZZ'))]
    dy_mc = mc[(mc.mc_cat == "DY")]

    # double weights, since we're only using half for each
    diboson_mc.weight *= 2
    dy_mc.weight *= 2

    # choose half of the diboson MC to use
    template_indices = numpy.random.choice([True, False], size=diboson_mc.weight.count(), replace=True)
    template_mc = diboson_mc[template_indices]

    # add the rest to the test sample
    test_mc = test_mc.append(diboson_mc[~template_indices])

    # now do the same for DY
    template_indices = numpy.random.choice([True, False], size=dy_mc.weight.count(), replace=True)
    template_mc = template_mc.append(dy_mc[template_indices])

    # add the rest to the test sample
    test_mc = test_mc.append(dy_mc[~template_indices])


    f = open("tables/table{}.tex".format(flavor), 'w')

    f.write("""
    \\begin{tabular}{|l|c|c|c|c|c|c|c|c|}
    """)

    table = defaultdict(lambda : defaultdict(lambda : defaultdict(dict)))
    for cut in numpy.arange(70., 130., 10.) :
        predictions = BackgroundFit.do_bkg_fit( test_mc, template_mc, cut, flavor, 0.07)
        truth = mctruth.get_mc_truth(test_mc, cut, flavor)

        for mctype in predictions['low'].keys() :
            table['Low \mctp\ Region'][mctype]['prediction'][cut] = predictions['low'][mctype]
            table['Low \mctp\ Region'][mctype]['truth'][cut] = truth['low'][mctype]


        for mctype in predictions['high'].keys() :
            table['High \mctp\ Region'][mctype]['prediction'][cut] = predictions['high'][mctype]
            table['High \mctp\ Region'][mctype]['truth'][cut] = truth['high'][mctype]

    mc_nicename = {"Top":"Top", "VV":"Diboson", "DY":"\PZz/\Pgg$^*$", "W":"Fake Lepton", "Total":"Total"}

    mctypes=['Top', 'VV', 'DY', 'W', "Total"]

    regions = table.keys()
    regions.reverse()
    for region in regions :
        f.write("\hline\n \multicolumn{7}{|c|}{"+region+"}\\\\ \n \hline\n")
        mccutstring = ""
        for mctype in mctypes :
            if mctype in table[region].keys():
                cuts = table[region][mctype]['prediction'].keys()
                cuts.sort()
                if mccutstring == "" :
                    mccutstring += "\mctp\ cut "
                    for cut in cuts :
                        mccutstring += (" & "+str(cut))
                    mccutstring += "\\\\ \n \hline\n"
                    f.write(mccutstring)
                truthline = mc_nicename[mctype]+" truth"
                for cut in cuts :
                    truthline += " & ${:.2f}$".format(table[region][mctype]['truth'][cut])
                truthline += "\\\\ \n\hline\n"
                f.write(truthline)

                predline = mc_nicename[mctype]+" prediction"
                for cut in cuts :
                    predline += " & ${:.2f} \pm {:.2f}$".format(*table[region][mctype]['prediction'][cut])
                predline += "\\\\ \n \hline\n"
                f.write(predline)

    f.write("""\\end{tabular}
    """)

    f.close()

if __name__ == '__main__':
    run_mc_closure('_of')
