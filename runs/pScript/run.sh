parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
echo "Plotting..."
#clear
if [ $# -eq 0 ]
  then
    time_start=0
    time_end=53
else
    time_start="$1"
    time_end="$2"
fi
dif="$(($time_end-$time_start))"

#reg
#python3 skeleTrace.py
#ffmpeg -r 72 -y -threads 4 -i 2DEvolve/skeleplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4

#Compare
#python3 vofsmooth.py
#ffmpeg -r 4 -y -threads 4 -i 2DEvolve/svofplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4

#inttrace
python3 splinecalc.py $time_start $time_end

if [ $dif -eq 0 ]
  then
    echo "cant animate because time error"
else
    ffmpeg -r 5 -y -threads 4 -i 2DEvolve/skele2intplt-%03d.png -pix_fmt yuv420p SplineEvolve.mp4
fi
