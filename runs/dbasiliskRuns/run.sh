parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
clear
./clean.sh
./compile.sh
./drop
./../pScript/run.sh
