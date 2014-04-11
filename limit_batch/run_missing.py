import json
import glob
import re
from collections import defaultdict
import envoy

with open("../slep_masses.json") as f:
    masses = json.load(f)

masses = [tuple(map(int, mass)) for mass in masses if mass[0] < 500 and mass[0] - mass[1] > 50]

def get_masses(filename):
    return tuple(map(int, re.search("_(\d+)_(\d+)[_\.]", filename).groups()))

def find_missing():
    filenames = glob.glob("slep_combined/*.root")
    finished_masses = map(get_masses, filenames)
    missing = set(masses) - set(finished_masses)
    return missing

missing = find_missing()

def find_missing_lines(missing):
    mass_to_line = defaultdict(list)
    missing_lines = []
    with open("slep_limits_new.txt") as f:
        for i, line in enumerate(f):
            masses = get_masses(line)
            mass_to_line[masses].append(i)
    for mass in missing:
        missing_lines.append(mass_to_line[mass])

    return missing_lines

missing_lines = find_missing_lines(missing)

def line_to_job(line):
    return (line-1)*20 + 1

for group in missing_lines:
    low = line_to_job(min(group))
    high = line_to_job(max(group)+20)
    command = "bsub -J 'slep_limits[{0}-{1}]' -q 1nh slep_limits_lsf.sh".format(low, high)
    envoy.run(command)