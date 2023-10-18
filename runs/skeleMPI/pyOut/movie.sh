ffmpeg -framerate 15 -pattern_type glob -i 'Img/plot-*.png' -c:v libx264 -pix_fmt yuv420p Movie/output.mp4
