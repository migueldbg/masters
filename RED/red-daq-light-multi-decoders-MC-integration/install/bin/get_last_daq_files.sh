#!/bin/bash

DATAPATH=/data/ReD-data/data/
TEMP_FILE_LIST=temp_last_daq_files.txt
EXEC=$1
#~/ddurso/red-daq-inst/bin/viewer_online
RUN=$2
if [ -z "$1" ] || [ -z "$2" ]
  then
    echo "No argument supplied"
fi

if [ -e $TEMP_FILE_LIST ]
then
    echo "removing temp list "
    rm $TEMP_FILE_LIST
fi
touch $TEMP_FILE_LIST

echo $DATAPATH
for board in `seq 0 3`;
do
    echo "board " $board
    if [ ! -z "$2" ] ; then 
	RUN=$2
	filename=` find $DATAPATH  -name "run_${RUN}_b0${board}*" |  xargs --no-run-if-empty ls --full-time  | awk '{ print $9 }' | sort | grep run | tail -1`
	echo " RUN " $RUN
    else
	filename=` find $DATAPATH  -name "run_*_b0${board}*" |  xargs --no-run-if-empty ls --full-time  | awk '{ print $9 }' | sort | grep run | tail -1` 
    fi
    echo $filename
    [ -z "$filename" ]  && continue; 
    echo "$board $filename " >> $TEMP_FILE_LIST

done

cat $TEMP_FILE_LIST

CMD="$EXEC -l  $TEMP_FILE_LIST "
$CMD

