#!/bin/bash

#for set in A B C D
for set in E F G H I J K L
do
  for mode in 0p0 0p1 0p2 0p3 0p4 0p5 0p6 0p7 0p8 0p9 1p0
  do
    ./make_script_2HDM.sh        $mode $set
    #./make_script_2HDM_double.sh $mode $set
  done
done
