import simplejson as json
from math import floor, log10

with open("fit_results.json") as f:
    results = json.load(f)

with open("mc_preds.json") as f:
    mc = json.load(f)

name_to_nicename = {'top': 'Top',
                    'vv': 'Diboson and Rare SM',
                    'z': 'Z/$\Pgg^*$',
                    'fake': 'Non-prompt',
                    'sum': 'Sum'}

of = results['of']
sf = results['sf']
of_mc = mc['of']
sf_mc = mc['sf']

def format_result(tpl):
    val, err = tpl
    digits = -int(floor(log10(err)))
    digits = max(digits, 0)
    if digits == 0:
        val, err = int(round(val, digits)), int(round(err, digits))
    if val > err:
        s = "${0}\pm${1}$".format(val, err)
    else:
        s = "${0}^{{+{1}}}_{{-{2}}}$".format(val, err, val)
    return s

for channel in ['top', 'vv', 'z', 'fake', 'sum']:
    # s = name_to_nicename[channel]
    s = name_to_nicename[channel]
    s += " & "
    s += format_result(of[channel])
    s += " & "
    s += str(int(round(of_mc[channel], 0)))
    s += " & "
    s += format_result(sf[channel])
    s += " & "
    s += str(int(round(sf_mc[channel], 0)))
    s += "\\\\"
    print s
