#! /bin/bash

source $HOME/.bash_profile

cd $HOME/work/MCTSusy/limit_batch

k5reauth -p neggert -k ~/user.keytab -- python -u get_all_limits.py chi_limits.txt chi_combined_20140128 8nh ../limits/chi_templates_*_combined_meas_model.root
