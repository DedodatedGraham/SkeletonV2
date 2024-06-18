#!/bin/sh

clear

max_level=13
rhog=5.00E-2
mul=4.84E-3
mug=4.84E-5
u0=127.65
#static
L=16
t_out=0.001
t_end=0.5
hx=-3.0
hy=0.0
femax=1e-3
uemax=1.225E+00
maxruntime=14400
tamr=0.


./clean.sh
./compile.sh
#./drop $max_level $L $u0 $t_out $t_end $hx $hy $rhog $mul $mug $femax $uemax $maxruntime $tamr > drop.out 2> drop.err
mpirun -np 12 ./drop $max_level $L $u0 $t_out $t_end $hx $hy $rhog $mul $mug $femax $uemax $maxruntime $tamr 2> drop.err
#mpirun -np 12 xterm -e gdb -ex run --args ./drop $max_level $L $u0 $t_out $t_end $hx $hy $rhog $mul $mug $femax $uemax $maxruntime $tamr 2> drop.err


echo
echo "Job finished at `date`"

