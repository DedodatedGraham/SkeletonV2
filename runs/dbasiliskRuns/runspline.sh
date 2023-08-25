start=`date +%s.%N`
tstart=0
tend=39
./clean.sh
./../pScript/clean.sh
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
clear
./compile.sh
if [ $# -eq 0 ]
  then
    echo "No Input"
    tstart=0
    tend=39
else
    tstart="$1"
    tend="$2"
fi
for D in $(seq 1 1 8)
do
    save="N$D"
    mdir="$save"
    saveForm="-s.mp4"
    save="$save$saveForm"
    mdir=$(echo $mdir | tr -d '\r')
    echo $mdir
    echo $save
    mkdir $mdir
    echo $(pwd)
    cd "$mdir"
    echo $(pwd)
    ./../drop $tstart $tend $D ../../basiliskRuns/dump-
    cd ../
    ./../pScript/runspline.sh $save 0 $mdir $tstart $tend &
done
end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
echo $runtime
