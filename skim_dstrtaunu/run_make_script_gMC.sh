#!/bin/bash

for stream in 0 1 2 3 4 5
#for stream in 6 7 8 9 # valid only for mixed and charged
do
    bsub -q e ./make_script_gMC.sh      uds     $stream
    bsub -q e ./make_script_gMC.sh      charm   $stream
    bsub -q e ./make_script_gMC.sh      mixed   $stream
    bsub -q e ./make_script_gMC.sh      charged $stream
done

