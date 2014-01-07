import simplejson as json
from math import floor, log10
from print_results_table import format_result

with open("fit_results2.json") as f:
    results = json.load(f)

with open("mc_preds2.json") as f:
    mc = json.load(f)

name_to_nicename = {'top': 'Flavor Symmetric',
                    'vv': 'Non-FS Diboson',
                    'z': 'Z/$\Pgg^*$',
                    'fake': 'Non-prompt',
                    'sum': 'Sum'}

sf = results['sf']
sf_mc = mc['sf']


def main():
    for channel in ['top', 'vv', 'z', 'fake', 'sum']:
        # s = name_to_nicename[channel]
        s = name_to_nicename[channel]
        s += " & "
        s += format_result(sf[channel])
        s += " & "
        s += str(int(round(sf_mc[channel], 0)))
        s += "\\\\"
        print s

if __name__ == '__main__':
    main()
