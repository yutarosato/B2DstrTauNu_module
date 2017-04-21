#!/bin/bash



#for set in A B C D E F G H I J K L M N O P
for set in AA BB CC DD EE FF GG HH
do
    bsub -q e ./make_script_sigMC.sh sig_kpi $set
    bsub -q e ./make_script_sigMC.sh nor_kpi $set
done

#bsub -q e ./make_script_sigMC.sh sig_kpi_dpi0 A
#bsub -q e ./make_script_sigMC.sh nor_kpi_dpi0 A
#bsub -q e ./make_script_sigMC.sh sig_kpipi0   A
#bsub -q e ./make_script_sigMC.sh nor_kpipi0   A

#bsub -q e ./make_script_sigMC.sh nor_kpi_dstpipi AA
#bsub -q e ./make_script_sigMC.sh nor_kpi_dstpipi BB
#bsub -q e ./make_script_sigMC.sh nor_kpi_dstst   AA
#bsub -q e ./make_script_sigMC.sh nor_kpi_dstst   BB
