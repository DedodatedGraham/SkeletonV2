#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"


rm -r *.dat 2> /dev/null
rm -r dumps/dump-* 2> /dev/null
rm -r dumps/snapshot-* 2> /dev/null
rm -r .qcc* 2> /dev/null
rm -r Img/*.png 2> /dev/null
rm -r drop 2> /dev/null
