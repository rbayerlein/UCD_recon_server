#!/bin/bash

# this is preliminary script to monitor the folder size of a given folder for 1000sec in steps of delta_t
PWD=$(pwd)
echo "this dir: $PWD"

folder=$1
cd $folder 
folder_absPath=$(pwd)

#place log file 1 level higher
cd ..
dest_dir=$(pwd)
cd $PWD

LOGFILENAME=$dest_dir/"monitorFolder.log"
echo " LOGFILENAME: $LOGFILENAME"

if test -f $LOGFILENAME; then
	echo "overwriting existing log file"
	rm $LOGFILENAME
fi

delta_t=10	# time increment in between steps

## begin main program

for (( i = 0; i < 1000; i++ )); do
	CHECK=$(du -sb $folder_absPath | cut -f1)
	let curr_t=$i*$delta_t
	echo -e "$curr_t \t $CHECK" >> $LOGFILENAME
	sleep $delta_t
done
