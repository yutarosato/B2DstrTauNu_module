#! /bin/tcsh -f
if( $#argv != 4 )then
    echo " Usage : $0 [que-type] [nor/sig] [mode] [set]"
    echo '[mode] = kpi kpipi0 kpi_dpi0'
    exit 1
endif


bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_07_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_09_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_11_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_13_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_15_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_17_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_19a_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_19b_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_21_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_23_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_25_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_27_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_31_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_33_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_35_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_37a_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_37b_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_39_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_41a_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_41b_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_43_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_45a_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_45b_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_47_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_49_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_51_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_55_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_61_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_63_set${4}.sh
bsub -q ${1} script/${2}_${3}/sigMC_${2}_${3}_65_set${4}.sh
