start=`date +%s.%N`
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
clear
./clean.sh
./../pScript/clean.sh
./compile.sh
#valgrind ./drop 0 39 3
./drop 0 39 3
./../pScript/run.sh
end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
echo $runtime
