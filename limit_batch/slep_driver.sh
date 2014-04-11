#! /bin/bash

source $HOME/.bash_profile

cd $HOME/work/MCTSusy/limit_batch

k5reauth -p neggert -k ~/user.keytab -- python -u get_all_limits.py slep_limits.txt slep_combined 1nh ../limits/slep_templates_*_combined_meas_model.root
