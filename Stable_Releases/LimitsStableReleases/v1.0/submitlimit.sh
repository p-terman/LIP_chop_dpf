#!/bin/bash
#PBS -q regular
#PBS -l nodes=5:ppn=4
#PBS -l walltime=20:00:00
#PBS -N limit_2_30phe
#PBS -e limit_2_30phe.$PBS_JOBID.err
#PBS -o limit_2_30phe.$PBS_JOBID.out
#PBS -m ae
 
cd $PBS_O_WORKDIR
#pbsdsh $PBS_O_WORKDIR/runlimit.sh $PBS_VNODENUM
./runlimit.sh $PBS_VNODENUM
