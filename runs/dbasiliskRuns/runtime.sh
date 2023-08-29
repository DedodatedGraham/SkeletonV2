start=`date +%s.%N`
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
clear
if [ $# -eq 0 ]
  then
    time_start=0
    time_end=39
    delta=3
else
    if [ $# -eq 1 ]
      then
        time_start="$1"
        time_end=39
        delta=3
    else
        if [ $# -eq 2 ]
          then
            time_start="$1"
            time_end="$2"
            delta=3
        else
            time_start="$1"
            time_end="$2"
            delta="$3"
        fi
    fi
fi
./clean.sh
./../pScript/clean.sh
./compile.sh
#valgrind ./drop $time_start $time_end $delta
./drop $time_start $time_end $delta
./../pScript/run.sh $time_start $time_end
end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
echo $runtime
