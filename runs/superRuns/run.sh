#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -o ds_job.out
#PBS -e ds_job.err
#PBS -q ling2

module load openmpi

echo "------------------"
echo
echo "Job working directory: $PBS_O_WORKDIR"
echo

num=`cat $PBS_NODEFILE | wc -l`
echo "Total processes: $num"
echo

echo "Job starting at `date`"
echo

cd $PBS_O_WORKDIR

#serial
mpiexec -n $num drop 1>tmpout 2>log
mpiexec -n 1 -machinefile $PBS_NODEFILE python3 ../pScript/skeleTrace.py

#openMPI
#mpiexec -n $num drop_mpi 1>tmpout 2>log

echo
echo "Job finished at `date`"
