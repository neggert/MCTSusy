import json

f = open("Cert_190456-201678_8TeV_PromptReco_Collisions12_JSON_remove2012C-v1_197770_198913.txt")

rl = json.load(f)

my_rl = {}

for key in rl.keys():
    runs = range(rl[key][0][0], rl[key][0][1]+1)
    for r in rl[key][1:]:
        runs.extend(range(r[0], r[1]+1))
    my_rl[int(key)] = runs


def good_lumi(event):
    try :
        return event['lumi'] in my_rl[event['run']]
    except KeyError:
        return False

def good_lumi_old(event):
    try :
        return check_range(event['lumi'], rl[str(event['run'])])
    except KeyError:
        return False
