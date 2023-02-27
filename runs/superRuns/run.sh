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

#openMPI
#mpiexec -n $num drop_mpi 1>tmpout 2>log

#animate doesnt work right now from here?
#mpiexec -n 1 -machinefile $PBS_NODEFILE python3 ../pScript/skeleTrace.py
#mpiexec -n 1 -machinefile $PBS_NODEFILE ffmpeg -r 72 -y -threads 4 -i ../pScript/2DEvolve/skeleplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4

echo
echo "Job finished at `date`"
