parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
#clear
if [ $# -eq 0 ]
  then
    save_name="2DEvolve.mp4"
else
    save_name="$1"
    if [ $2 -eq 1 ]
        then
            ./clean.sh
    fi
fi
echo $save_name
#reg
#python3 skeleTrace.py
#ffmpeg -r 72 -y -threads 4 -i 2DEvolve/skeleplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4

#Compare
#python3 vofsmooth.py
#ffmpeg -r 4 -y -threads 4 -i 2DEvolve/svofplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4

#inttrace
python3 splinecalc.py


ffmpeg -r 5 -y -threads 4 -i 2DEvolve/skele2intplt-%03d.png -pix_fmt yuv420p $save_name &
