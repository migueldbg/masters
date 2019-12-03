#!/bin/bash

DATAPATH=/NAS/data_from_daq1/rawdata/
COUNTER=$1
SER=$3
while [  $COUNTER -le $2 ]; do
#  echo item: $COUNTER
  ./makelist_multi_naples.sh $COUNTER
  mv run_$COUNTER.lst ../list  
  ../RedLevel1 -o $HOME/rootfiles/run_$COUNTER.root -l ../list/run_$COUNTER.lst -c ../cfg/global_ph2.cfg -m ../cfg/channelmapping_Na.cfg -e ../ser/ser_$SER.cfg
  let COUNTER=COUNTER+1
done
