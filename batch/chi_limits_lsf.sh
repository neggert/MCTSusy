#! /bin/sh
#BSUB -J "neggertArray[1-521]"
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -q 1nw

echo $LSB_JOBINDEX



cd $HOME/work/CMSSW_5_2_5/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.04/x86_64-slc5-gcc46-opt/root/
source bin/thisroot.sh
pythonbrew use 2.7.3


cd $HOME/work/MCTSusy
pythonbrew py -p 2.7.3 set_limits.py batch sig_chi.root chi_masses.json $(($LSB_JOBINDEX-1)) limits/results/chi_$LSB_JOBINDEX.dat --ncpu=8
