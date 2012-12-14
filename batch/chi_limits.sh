#PBS -t 0-130

#!/bin/sh

source /home/uscms33/setup.sh
source /home/uscms33/setupROOT.sh

cd /home/uscms33/MCTSusy/

echo $PBS_ARRAYID

./set_limits.py batch sig_chi.root chi_masses.json $PBS_ARRAYID limits/results/slep_$PBS_ARRAYID.dat
