#!/bin/bash
export ROOTSYS=/homes/lmc/Run3Limits/root
export PATH=${ROOTSYS}/bin:${PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib
source $ROOTSYS/bin/thisroot.sh

#mass=( 6 7 8 9 10 12 14 17 21 27 33 40 50 60 70 80 100 300 1000 2000 )
#echo -n "Running mass = "
#echo ${mass[$1]}

#wdir=$(printf "mass%i_output" ${mass[$1]})     # dir00, dir01, dir02,...
wdir=$(printf "mass%i_output" $1)     # dir00, dir01, dir02,...
mkdir $wdir

#./runMass wsData_20131002.root ${mass[$1]}
./runMass wsData_20131002.root $1
mv *.txt $wdir
mv ws_model_*.root $wdir
mv *.pdf $wdir

