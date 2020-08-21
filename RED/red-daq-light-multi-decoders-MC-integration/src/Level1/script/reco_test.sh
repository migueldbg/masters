#!/bin/bash
### Example two phases : reco.sh 528 ph2 513
###        single phase: reco.sh 516 ph1 499
###        laser run:    reco.sh 499 laser
### Example to run reco on two phase run from 529 to 535:
###                      for _run in `seq 529 535 `; do echo $_run ; ./reco.sh $_run ph2 513  ; done  
### WARNING: add a new env variable if you want to use a different channel mapping file!!!!!!
### Simone 10-07-18

runumber=$1
globalnumber=$2
sernumber=$3
ch_mapping=$4
list_name="run_$runumber.lst"
global_name="global_$globalnumber"
laser_name="laser"
ser_name="ser_$sernumber"
RECODIR=/storage/DATA-02/darkside/red/reco/rm3reco/lns
file_name="run_${runumber}_${laser_name}"
file_name1="run_${runumber}"
file_name2="run_${runumber}_${global_name}_${ser_name}"
# Replace here, in order to process on the LNS DAQ server
#SCRIPTDIR=/home/pandola/red-daq-light/src/Level1/script
LEVEL1DIR=/storage/local/home/darkside/nrossi/Level1
SCRIPTDIR=/storage/local/home/darkside/nrossi/Level1/script 

cd $LEVEL1DIR
#./script/makelist_multi_lns.sh $runumber
if [ ${globalnumber:0:5} == laser ]
	then
	./RedLevel1 -l $list_name -m cfg/$ch_mapping.cfg -c cfg/$globalnumber.cfg -o $RECODIR/$file_name.root --laser
	echo " Reco $file_name done! "
else
   ./RedLevel1 -l $list_name -m cfg/$ch_mapping.cfg -c cfg/$global_name.cfg -e ser/$ser_name.cfg -o $RECODIR/$file_name1_$2.root 
echo " Reco $file_name1 done! ";
fi

cd $SCRIPTDIR
