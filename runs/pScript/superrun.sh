#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -o fig.out
#PBS -e fig.err

cd $PBS_O_WORKDIR
echo "Starting at $PBS_O_WORKDIR"
echo "node file $PBS_NODEFILE"
num=`cat $PBS_NODEFILE | wc -l`

#Skeleton
#mpiexec -n 1 -machinefile $PBS_NODEFILE python3 skeleTrace.py
#mpiexec -n 1 -machinefile $PBS_NODEFILE ffmpeg -r 10 -y -threads 4 -i 2DEvolve/skeleplt-%03d.png -pix_fmt yuv420p -vf scale=-1:720 2DEvolve.mp4

#Normal
mpiexec -n 1 -machinefile $PBS_NODEFILE python3 normtest.py
mpiexec -n 1 -machinefile $PBS_NODEFILE ffmpeg -r 10 -y -threads 4 -i 2DEvolve/normplt-%03d.png -pix_fmt yuv420p -vf scale=-1:720 2DEvolve.mp4

echo
echo "Job finished at `date`"
