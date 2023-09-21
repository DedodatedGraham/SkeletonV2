./clean.sh
./compile-ns.sh
export OMP_NUM_THREADS=28
#valgrind ./_drop 10
./_drop 10
rm _drop
