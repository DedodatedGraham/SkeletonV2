start=`date +%s.%N`
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
./compile.sh
for L in $(seq 5 1 9)
do
    echo $L
    for D in $(seq 1 2 9)
    do
        ./cleanspline.sh
        echo $D
        #clear
        ./drop $L $D
        save="N3"
        saveL="L"
        saveL="$saveL$L"
        saveD="D"
        saveD="$saveD$D"
        saveForm=".mp4"
        save="$save$saveL$saveD$saveForm"
        echo save
        ./../pScript/runspline.sh $save 0
    done
done
end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
echo $runtime
