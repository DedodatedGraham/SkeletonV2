parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
./clean.sh
#reg
#python3 skeleTrace.py
#ffmpeg -r 72 -y -threads 4 -i 2DEvolve/skeleplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4

#Compare
python3 vofsmooth.py
ffmpeg -r 4 -y -threads 4 -i 2DEvolve/svofplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4
