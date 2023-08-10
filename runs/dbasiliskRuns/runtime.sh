start=`date +%s.%N`
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
clear
if [ $# -eq 0 ]
  then
    time_start=0
    time_end=39
else
    if [ $# -eq 1 ]
      then
        time_start="$1"
        time_end=39
    else
        time_start="$1"
        time_end="$2"

    fi
fi
./clean.sh
./compile.sh
#valgrind ./drop $time_start $time_end
./drop $time_start $time_end
#./../pScript/run.sh $time_start $time_end
end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
echo $runtime
