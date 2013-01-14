"""
Utility for filtering events based on good run list
"""

import json

f = open("Merged_190456-208686_8TeV_PromptReReco_Collisions12.json")

rl = json.load(f)

my_rl = {}

for key in rl.keys():
    runs = range(rl[key][0][0], rl[key][0][1]+1)
    for r in rl[key][1:]:
        runs.extend(range(r[0], r[1]+1))
    my_rl[int(key)] = runs


def good_lumi(event):
    """Return true if event has a good run and lumi, false otherwise"""
    try :
        return event['lumi'] in my_rl[event['run']]
    except KeyError:
        return False

def good_lumi_old(event):
    try :
        return check_range(event['lumi'], rl[str(event['run'])])
    except KeyError:
        return False
