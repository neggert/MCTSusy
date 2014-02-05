#! /bin/sh

FILE=$1
VAL=$2

set --


source /afs/cern.ch/sw/lcg/external/gcc/4.6.3/x86_64-slc6/setup.sh

source $HOME/mctpy/bin/activate

source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.04/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh

printenv

echo $(which python)

kinit neggert@CERN.CH -k -t $HOME/user.keytab

cd $HOME/work/MCTSusy/limit_batch
python run_one.py $FILE $VAL $LSB_JOBID $HOME/work/MCTSusy/limit_batch/tmp/$LSB_JOBID.root

rm core.*
