#! /bin/sh
#BSUB -J "chi_limits[1-500]`
#BSUB -q 8nh

source /afs/cern.ch/sw/lcg/external/gcc/4.6.3/x86_64-slc6/setup.sh

source $HOME/mctpy/bin/activate

source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.04/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh

cd $HOME/work/MCTSusy/limit_batch
python run_one.py batch chi_limits_extra10.txt $(($LSB_JOBINDEX-1)) chi chi_limits_extra
