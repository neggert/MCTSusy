#! /bin/sh

for i in {1..5};
do
	bsub -J "chi_limits[$(((i-1)*1000+1))-$((i*1000))]" < chi_limits_lsf.sh
done