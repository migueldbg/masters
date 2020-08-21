#! /bin/bash

grep Evt /data/ReD-data/data/run_$1/log/run_$1_b0*.log | head > test.txt
awk '{if (NR!=1&&NR<12) t0[NR-2]=$5; else if (NR!=12&&NR<23) t1[NR-13]=$5; else if (NR!=23&&NR<34) t2[NR-24]=$5}; END {for (j=0;j<10;j++) print t0[j], t1[j], t2[j]}'  test.txt
