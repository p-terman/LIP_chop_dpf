#!/bin/bash

#This script is to find the paths of the evt and  100 pulse rq files for the
#select reprocessing dpf for LIP pulse chopping. And then passes (and runs) that
#information to the python script (which finds the LUG info and runs the Matlab script)
# arg 1 is the form luxAA_yyyymmddThhmm
# arg 2 LUXcode_path i.e. the path the the modified modules LUXcode
# arg 3 is data_processing_xml_path
# note that nothing is really being done with arg2 and 3. They are just along for the ride. 
# This program only works assuming that the 100 pulse rq are the highest cp version of an rq  
#
# version 1.0 20190220 PAT

evtPath=$($(./get_dataset_paths_alex.sh $1 evt) 
rqList=$(./get_dataset_paths_alex.sh $1 rq) #make sure you know where you are first
stringArray=($rqList)

a=${#stringArray[@]}
for f in $(seq 0 1 $((10#$a - 1))) #this is the length of the array

    do

        numArray[$f]=${stringArray[$f]##*_cp}
        numArray[$f]=$((10#${numArray[$f]}))

done

max=${numArray[0]}
max_ind=$((0))
for i in $(seq 1 1 $((${#numArray[@]} - 1)))
    do

    if test ${numArray[$i]} -gt $max
    then
        max=${numArray[$i]}
        max_ind=$((10#$i))
    fi

done

rqPath= ${stringArray[$max_ind]}

python run_select_dpf_brown.py $1 $evtPath $rqPath $2 $3


