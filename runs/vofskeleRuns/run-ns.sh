#Compile first case of drop
./compile-ns.sh
export OMP_NUM_THREADS=28
./_drop 11
rm _drop
