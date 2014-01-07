import simplejson as json
from math import floor, log10
from print_results_table import format_result

with open("count_fit_results.json") as f:
    results = json.load(f)

name_to_nicename = {'top': 'Top',
                    'vv': 'Diboson and Rare SM',
                    'z': 'Z/$\Pgg^*$',
                    'fake': 'Non-prompt',
                    'sum': 'Sum'}

of = results['of']
sf = results['sf']


for channel in ['top', 'vv', 'z', 'fake', 'sum']:
    # s = name_to_nicename[channel]
    s = name_to_nicename[channel]
    s += " & "
    s += format_result(of[channel]['low'])
    s += " & "
    s += format_result(of[channel]['high'])
    s += " & "
    s += format_result(sf[channel]['low'])
    s += " & "
    s += format_result(sf[channel]['high'])
    s += "\\\\"
    print s
