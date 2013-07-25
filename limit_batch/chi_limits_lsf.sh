#! /bin/sh
#BSUB -J "chi_limits[1-500]`
#BSUB -q 8nh


echo $LSB_JOBINDEX

cd $HOME/work/CMSSW_5_2_5/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.04/x86_64-slc5-gcc46-opt/root/
source bin/thisroot.sh
pythonbrew use 2.7.3


cd $HOME/work/MCTSusy/limit_batch
pythonbrew py -p 2.7.3 run_one.py batch chi_batch_extra8.txt $(($LSB_JOBINDEX-1)) chi extra8
