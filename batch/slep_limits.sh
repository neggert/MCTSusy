#PBS -t 0-71

#!/bin/sh

source /home/uscms33/setup.sh
source /home/uscms33/setupROOT.sh

cd /home/uscms33/MCTSusy/

echo $PBS_ARRAYID

./set_limits2.py batch signal_slep.root slep_masses.json $PBS_ARRAYID limits/results/slep_$PBS_ARRAYID.dat --channels=sf
