#ffmpeg -framerate 30 -pattern_type glob -i 'Img/plot-*.png' -c:v libx264 -pix_fmt yuv420p Movie/output.mp4


ffmpeg -framerate 15 -pattern_type glob -i 'Img/plot-*.png' -c:v libx264 -pix_fmt yuv420p Movie/skeletonthin$1.mp4


