parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
rm -r 2DEvolve/skeleplt*.png 2> /dev/null
rm -r 2DEvolve/skele2intplt*.png 2> /dev/null
rm -r 2DEvolve/splinecalcplt*.png 2> /dev/null
rm -r 2DEvolve/normplt*.png 2> /dev/null
rm -r 2DEvolve/compplt*.png 2> /dev/null
rm -r 2DEvolve/svofplt*.png 2> /dev/null
rm *.mp4 2> /dev/null
rm fig.* 2> /dev/null
rm -r N*D* 2> /dev/null
