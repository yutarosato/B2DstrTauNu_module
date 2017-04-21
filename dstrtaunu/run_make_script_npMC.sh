#!/bin/bash

for PAR in 0p0
do
  #for STREAM in $(seq 0 39)
  for STREAM in $(seq 10 39)
  do
    ./make_script_npMC.sh ${PAR} ${STREAM}
  done
done

exit
  
for PAR in 0p5
do
  for STREAM in $(seq 0 9)
  do
    ./make_script_npMC.sh ${PAR} ${STREAM}
  done
done
