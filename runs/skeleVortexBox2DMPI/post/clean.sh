#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"


rm -r *.dat 2> /dev/null
rm -r .qcc* 2> /dev/null
rm -r Img/*.png 2> /dev/null
rm -r Movie/*.mp4 2> /dev/null
rm -r drop 2> /dev/null
rm -r dat/* 2> /dev/null
rm -r mtrace* 2> /dev/null
rm -r core 2> /dev/null
