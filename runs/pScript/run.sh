./clean.sh
python3 skeleTrace.py
ffmpeg -r 72 -y -threads 4 -i 2DEvolve/skeleplt-%03d.png -pix_fmt yuv420p 2DEvolve.mp4
