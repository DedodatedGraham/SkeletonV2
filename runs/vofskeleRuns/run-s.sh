./clean.sh
./compile-s.sh
export OMP_NUM_THREADS=28
#valgrind ./_drop 10 335
./_drop 10

#run in mpi
#mpirun -np 2 _drop

