#! /bin/sh

for i in {1..3};
do
	bsub -J "chi_limits[$(((i-1)*1000+1))-$((i*1000))]" -q 8nh -o /dev/null < chi_limits_lsf.sh
	# bkill -J "chi_limits[$(((i-1)*1000+1))-$((i*1000))]"
done
