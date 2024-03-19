#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"


rm -r *.dat 2> /dev/null
rm -r dumps/dump-* 2> /dev/null
rm -r dumps/snapshot-* 2> /dev/null
rm -r .qcc* 2> /dev/null
rm -r Img/* 2> /dev/null
rm -r Movie/*.mp4 2> /dev/null
rm -r cyl 2> /dev/null
rm -r dat/* 2> /dev/null
rm  log* 2> /dev/null
rm -r mtrace* 2> /dev/null
rm -r valgrind.txt 2> /dev/null
rm -r core 2> /dev/null
rm -r cp-* 2> /dev/null
rm omega 2> /dev/null
rm perfs 2> /dev/null
