#!/bin/bash

#for set in A B
for set in C D E F G H I J K L
do
for type in R2 S1
do
#  for mode in 0p000 0p090 0p180 0p260 0p340 0p400 m0p040 m0p100 m0p150
#  for mode in 0p250 0p360
  for mode in m0p030
  do
    ./make_script_LQMC.sh        $mode $type $set
#    ./make_script_LQMC_double.sh $mode $type $set
  done
done
done
