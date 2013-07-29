import sys
from pandas import *


files = sys.argv[1:]

hdf_files = [f for f in files]

stores = []

for f in hdf_files:
    stores.append(HDFStore(f))

# combine all the dataframes
df = stores[0]['data']
for s in stores[1:]:
	try:
	    df = df.append(s['data'], ignore_index=True)
	except KeyError:
		raise


# drop duplicates
df.drop_duplicates(['run','lumi','event'], inplace=True)

import json

f = open("import/Merged_190456-208686_8TeV_PromptReReco_Collisions12.json")

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

good_events = df.apply(good_lumi, 1)

df = df[good_events]

import IPython
IPython.embed()