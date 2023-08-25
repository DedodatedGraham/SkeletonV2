parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
#clear
if [ $# -eq 0 ]
  then
    save_name="2DEvolve.mp4"
else
    save_name="$1"
    adjust_path="$3"
    if [ $2 -eq 1 ]
        then
            ./clean.sh
    fi
    time_start="$4"
    time_end="$5"
fi
echo $save_name
adjust_path=$(echo $adjust_path | tr -d '\r')
#reg
#python3 skeleTrace.py
#ffmpeg -r 72 -y -threads 4 -i 2DEvolve/skeleplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4

#Compare
#python3 vofsmooth.py
#ffmpeg -r 4 -y -threads 4 -i 2DEvolve/svofplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4

#inttrace


dif="$(($time_end-$time_start))"
if [ $dif -eq 0 ]
  then
    python3 splinecalc.py $time_start $time_end $adjust_path 1
else
    mkdir "$adjust_path"
    python3 splinecalc.py $time_start $time_end $adjust_path
    ffmpeg -r 5 -y -threads 4 -i $adjust_path/$adjust_path-%03d.png -pix_fmt yuv420p $save_name &
fi
