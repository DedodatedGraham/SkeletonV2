if [ $# -eq 0 ]
  then
    time_start=0
else
    time_start=$1
fi


./clean.sh
./compile-s.sh
export OMP_NUM_THREADS=28
#valgrind ./_drop 10
#./_drop 10 $time_start

#run MPI
mpirun -np 2 _drop 10 $time_start

rm _drop
