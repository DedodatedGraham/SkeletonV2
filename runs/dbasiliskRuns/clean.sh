parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
rm -r dump* 2> /dev/null
rm -r snapshot* 2> /dev/null
rm -r ../src 2> /dev/null
rm -r vof* 2> /dev/null
rm -r skeleton-* 2> /dev/null
rm -r intdata-* 2> /dev/null
rm -r smoothskeleton-* 2> /dev/null
rm -r reducedskeleton-* 2> /dev/null
rm -r nodePoint-* 2> /dev/null
rm -r boxDat-* 2> /dev/null
rm -r nodeDat-* 2> /dev/null
rm -r splinecalcDat-* 2> /dev/null
rm -r connectionDat-* 2> /dev/null
rm -r connectionidDat-* 2> /dev/null
rm -r splineBranchDat-* 2> /dev/null
rm -r angDat-* 2> /dev/null
rm -r .qcc* 2> /dev/null
rm -r infc* 2> /dev/null
rm -r mtrace* 2> /dev/null
rm drop 2> /dev/null
rm -r N*D*/ 2> /dev/null
