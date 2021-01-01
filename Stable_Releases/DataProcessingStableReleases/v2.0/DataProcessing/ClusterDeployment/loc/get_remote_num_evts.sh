#!/bin/bash

if [[ $# -gt 2 ]]; then
  echo "Error: Need dataset name. Add a -v for info"
  exit
fi


#NUM_EVTS=`ssh lux@gsk-41.het.brown.edu 'bash /Network/Servers/gsk-10.het.brown.edu/Volumes/gsk10sys/Users/lux/scripts/count_evt_files.sh $1' `
COMMAND="bash /Users/lux_mirror_user/scripts/count_evt_files.sh $1"
#echo COMMAND
#echo $COMMAND
NUM_EVTS=`ssh lux_mirror_user@gsk-41.het.brown.edu "$COMMAND" 2>/dev/null`
#echo NUM_EVENTS
#echo $NUM_EVTS >&2

if [[ $# == 2 ]]; then
  echo "$1 has $NUM_EVTS evt files"
else
  echo $NUM_EVTS
fi

