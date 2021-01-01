#!/bin/bash
#PBS -q debug
#PBS -l mppwidth=20
#PBS -l walltime=00:30:00
#PBS -N limit_2_30phe
#PBS -e limit_2_30phe.$PBS_JOBID.err
#PBS -o limit_2_30phe.$PBS_JOBID.out
#PBS -m ae
 
cd $PBS_O_WORKDIR
./runlimit.sh $PBS_VNODENUM
