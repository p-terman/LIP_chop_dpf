#!/bin/bash
#
# This script writes all dp started flags in every dataset in your evt
# directory.
#
# 2013-01-25 - JRV - Created
#

EVT_DIR="$HOME/scratch/evt/"
DP_STARTED_FLAG='ccv_dp_started'

dirs=`ls -1 $EVT_DIR`

for d in $dirs
do
    if [ -d $EVT_DIR/$d ]; then
        touch $EVT_DIR/$d/$DP_STARTED_FLAG
    fi  
done
