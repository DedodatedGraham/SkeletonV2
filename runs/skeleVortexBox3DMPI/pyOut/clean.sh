#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"


rm -r Img/*.png 2> /dev/null
rm -r Movie/*.mp4 2> /dev/null
rm -r __pycache__/ 2> /dev/null
