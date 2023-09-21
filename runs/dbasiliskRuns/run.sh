start=`date +%s.%N`
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
clear
./clean.sh
./../pScript/clean.sh
./compile.sh
valgrind ./drop 0 64 5
#./drop 0 53 5
#gdb -q -ex=r drop 0 50 8
./../pScript/run.sh
end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
echo $runtime
