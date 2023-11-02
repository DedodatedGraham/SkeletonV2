ffmpeg -framerate 30 -pattern_type glob -i 'Img/levelint-*.png' -c:v libx264 -pix_fmt yuv420p Movie/levelvof.mp4
ffmpeg -framerate 30 -pattern_type glob -i 'Img/ux-*.png' -c:v libx264 -pix_fmt yuv420p Movie/uxx.mp4
ffmpeg -framerate 30 -pattern_type glob -i 'Img/uy-*.png' -c:v libx264 -pix_fmt yuv420p Movie/uy.mp4





