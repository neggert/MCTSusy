import json

f = open("results_of.json")
results_of = json.load(f)
f.close()

f = open("results_sf.json")
results_sf = json.load(f)
f.close()

f = open("tables/results.tex", 'w')

f.write("""\\begin{table}
    \\begin{center}
    \caption{Results from shape fit to data in the low \mctp\ control region and extrapolation to the high \mctp\ signal region in both channels. Since the contributions of the individual backgrounds in the low \mctp\ region are derived from a fit to data, they by definition add to the number of events in that region in data.}
    \label{tab:fit}
    \\begin{tabular}{l|c|c|c|c}
        & \multicolumn{2}{c|}{Opposite Flavor} & \multicolumn{2}{c}{Same Flavor}\\\\
         & \\textbf{$5 < \mctp < 100 \GeV$} & \\textbf{$\mctp > 100 \GeV$} & \\textbf{$5 < \mctp < 100 \GeV$} & \\textbf{$\mctp > 100 \GeV$}\\\\
         \hline\n""")
f.write("Top & ${:.0f}\pm{:.0f}$ & ${:.1f}\pm{:.1f}$& ${:.0f}\pm{:.0f}$& ${:.1f}\pm{:.1f}$\\\\\n".format(results_of['low']['Top'][0], results_of['low']['Top'][1],
                                                                       results_of['high']['Top'][0], results_of['high']['Top'][1],
                                                                       results_sf['low']['Top'][0], results_sf['low']['Top'][1],
                                                                       results_sf['high']['Top'][0], results_sf['high']['Top'][1],
                                                                      ))
f.write("WW &  ${:.0f}\pm{:.0f}$ & ${:.1f}\pm{:.1f}$& ${:.0f}\pm{:.0f}$& ${:.1f}\pm{:.1f}$\\\\\n".format(results_of['low']['WW'][0], results_of['low']['WW'][1],
                                                                       results_of['high']['WW'][0], results_of['high']['WW'][1],
                                                                       results_sf['low']['WW'][0], results_sf['low']['WW'][1],
                                                                       results_sf['high']['WW'][0], results_sf['high']['WW'][1],
                                                                      ))
f.write("WZ &  ${:.0f}\pm{:.0f}$ & ${:.1f}\pm{:.1f}$& ${:.0f}\pm{:.0f}$& ${:.1f}\pm{:.1f}$\\\\\n".format(results_of['low']['WZ'][0], results_of['low']['WZ'][1],
                                                                       results_of['high']['WZ'][0], results_of['high']['WZ'][1],
                                                                       results_sf['low']['WZ'][0], results_sf['low']['WZ'][1],
                                                                       results_sf['high']['WZ'][0], results_sf['high']['WZ'][1],
                                                                      ))
f.write("ZZ &  ${:.1f}\pm{:.1f}$ & ${:.1f}\pm{:.1f}$& ${:.0f}\pm{:.0f}$& ${:.1f}\pm{:.1f}$\\\\\n".format(results_of['low']['ZZ'][0], results_of['low']['ZZ'][1],
                                                                       results_of['high']['ZZ'][0], results_of['high']['ZZ'][1],
                                                                       results_sf['low']['ZZ'][0], results_sf['low']['ZZ'][1],
                                                                       results_sf['high']['ZZ'][0], results_sf['high']['ZZ'][1],
                                                                      ))
f.write("DY &  ${:.0f}\pm{:.0f}$ & ${:.1f}\pm{:.1f}$& ${:.0f}\pm{:.0f}$& ${:.1f}\pm{:.1f}$\\\\\n".format(results_of['low']['DY'][0], results_of['low']['DY'][1],
                                                                       results_of['high']['DY'][0], results_of['high']['DY'][1],
                                                                       results_sf['low']['DY'][0], results_sf['low']['DY'][1],
                                                                       results_sf['high']['DY'][0], results_sf['high']['DY'][1],
                                                                      ))
f.write("W &   ${:.0f}\pm{:.0f}$ & ${:.1f}\pm{:.1f}$& ${:.0f}\pm{:.0f}$& ${:.1f}\pm{:.1f}$\\\\\n".format(results_of['low']['W'][0], results_of['low']['W'][1],
                                                                       abs(results_of['high']['W'][0]), results_of['high']['W'][1],
                                                                       results_sf['low']['W'][0], results_sf['low']['W'][1],
                                                                       abs(results_sf['high']['W'][0]), results_sf['high']['W'][1],
                                                                      ))
f.write("\hline\n")

f.write("Total & see caption & ${:.1f}\pm{:.1f}$ & see caption & ${:.1f}\pm{:.1f}$\\\\\n".format(results_of['high']['Total'][0], results_of['high']['Total'][1],
                                                                               results_sf['high']['Total'][0], results_sf['high']['Total'][1]
                                                                              ))
f.write("\hline\n\hline\n")
f.write("""Data & 2604 & 6 & 2438 & 6\\\\
    \end{tabular}
    \end{center}
\end{table}""")

f.close()

