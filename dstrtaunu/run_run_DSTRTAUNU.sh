#!/bin/bash

#for SET in A B C D E F G H I J K L M N O P
for SET in AA BB CC DD EE FF GG HH
#for SET in AA BB
do
  #./run_DSTRTAUNU_sigMC.sh s sig kpi $SET
  ./run_DSTRTAUNU_sigMC.sh s nor kpi $SET
  #./run_DSTRTAUNU_sigMC.sh s nor kpi_dstst   $SET
  #./run_DSTRTAUNU_sigMC.sh s nor kpi_dstpipi $SET
done

#./run_DSTRTAUNU_gMC_skim_all.sh 

#for STREAM in  0 1 2 3 4 5 # 6 7 8 9
#do
 #./run_DSTRTAUNU_gMC.sh      s uds     ${STREAM}
 #./run_DSTRTAUNU_gMC.sh      s charm   ${STREAM}
 #./run_DSTRTAUNU_gMC.sh      s mixed   ${STREAM}
 #./run_DSTRTAUNU_gMC.sh      s charged ${STREAM}
 #./run_DSTRTAUNU_gMC_skim.sh s uds     ${STREAM}
 #./run_DSTRTAUNU_gMC_skim.sh s charm   ${STREAM}
 #./run_DSTRTAUNU_gMC_skim.sh s mixed   ${STREAM}
 #./run_DSTRTAUNU_gMC_skim.sh s charged ${STREAM}
#done



