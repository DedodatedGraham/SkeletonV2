parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
./dbasiliskRuns/clean.sh
./pScript/clean.sh
