#! /bin/tcsh -f

if( $#argv != 2 )then
    echo " Usage : $0 [sig_kpi/nor_kpi/sig_kpipi0/nor_kpipi0/sig_kpi_dpi0/nor_kpi_dpi0] [set-name]"
    exit 1
endif

#****************************************************************************************************************
set MODE  = $1
set SET   = $2
set DIR   = '~/Store/dstrtaunu/mcprod/gsim/'$MODE'/mdst_set'$SET'/'
set MDL    = 'DSTRTAUNU'
set TAG   = '_set'$SET    # _NAME
set FLAG  = 'sigMC'
#****************************************************************************************************************
set RUN   = 'run_'$MDL'_'$FLAG'_'$MODE$TAG'.sh'
#****************************************************************************************************************

set list = `awk '{print $1}' ~/script/nBB.txt`

rm -f $RUN

# $1 : $DIR                                 # location of mdst file
# $2 : $SET($MDL/$TAG)                      # module name and etc...
# $3 : $f                                   # mdst file name for sigMC
# $4 : $DATA($FLAG)                         # MC flag
# $5 : $MODE($MODE)                         # module parameters
# $6 : $RUN                                 # run script file name

#  ************************ MC ***********************
echo '[ ' $MDL ' ' $FLAG ' ' ' ' $TAG ' ]' 

foreach f($list) 
  ./mk_script_sigMC.sh $DIR $MDL/$TAG $f $FLAG $MODE $RUN
end

chmod 744 $RUN
