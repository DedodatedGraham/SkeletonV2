ffmpeg -framerate 30 -pattern_type glob -i 'Img/levels-s-9-*.png' -c:v libx264 -pix_fmt yuv420p Movie/levels.mp4
ffmpeg -framerate 30 -pattern_type glob -i 'Img/sign-s-9-*.png' -c:v libx264 -pix_fmt yuv420p Movie/sign.mp4
ffmpeg -framerate 30 -pattern_type glob -i 'Img/svof-s-9-*.png' -c:v libx264 -pix_fmt yuv420p Movie/svof.mp4





