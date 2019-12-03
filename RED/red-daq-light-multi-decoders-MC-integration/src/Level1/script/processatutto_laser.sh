#!/bin/bash

DATAPATH=/NAS/data_from_daq1/rawdata/
COUNTER=$1
while [  $COUNTER -le $2 ]; do
#  echo item: $COUNTER
  ./makelist_multi_naples.sh $COUNTER
  mv run_$COUNTER.lst ../list  
  ../RedLevel1 -o $HOME/rootfiles/laser/run_$COUNTER.root -l ../list/run_$COUNTER.lst -c ../cfg/global_laser.cfg -m ../cfg/channelmapping_Na.cfg --laser -x
  let COUNTER=COUNTER+1
done
