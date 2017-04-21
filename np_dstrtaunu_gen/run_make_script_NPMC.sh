#!/bin/bash

#for set in A B
for set in C D E F G H I J K L
do


#for val in m5p00 m4p00 m3p00 m2p00 m1p00 0p00 1p00 2p00 3p00 # S1 : 9 points
#do
#    ./make_script_NPMC.sh         $val S1 $set
#    ./make_script_NPMC_double.sh  $val S1 $set
#done
  
#for val in m3p00 m2p00 m1p00  0p00  1p00 2p00 3p00 4p00 5p00 # S2 : 9 points
#do
#    ./make_script_NPMC.sh         $val S2 $set
#    ./make_script_NPMC_double.sh  $val S2 $set
#done

#for val in m2p50 m2p00 m1p50 m1p00 m0p50 0p00 0p50           # V1 : 7 points
for val in 0p10           # V1
do
    ./make_script_NPMC.sh         $val V1 $set
#    ./make_script_NPMC_double.sh  $val V1 $set
done

#for val in m0p50  0p00  0p50  1p00  1p50 2p00 2p50           # V2 : 7 points
for val in m0p10 # V2
do
    ./make_script_NPMC.sh         $val V2 $set
#    ./make_script_NPMC_double.sh  $val V2 $set
done

#for val in m0p15 m0p07  0p00  0p07  0p15 0p22 0p30 0p37 0p45 # T : 9 points
for val in m0p04 0p36 # T
do
    ./make_script_NPMC.sh         $val T $set
#    ./make_script_NPMC_double.sh  $val T $set
done





done





#    ./make_script_NPMC.sh         1p05 S1 $set
#    ./make_script_NPMC_double.sh  1p05 S1 $set
#    ./make_script_NPMC.sh        m2p10 V1 $set
#    ./make_script_NPMC_double.sh m2p10 V1 $set
#    ./make_script_NPMC.sh         1p88 V2 $set
#    ./make_script_NPMC_double.sh  1p88 V2 $set
