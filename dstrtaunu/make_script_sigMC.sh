#! /bin/tcsh -f

if( $#argv != 2 )then
    echo " Usage : $0 [sig_kpi/nor_kpi/sig_kpipi0/nor_kpipi0/sig_kpi_dpi0/nor_kpi_dpi0/nor_kpi_dstst/nor_kpi_dstpipi] [set-name]"
    exit 1
endif

#****************************************************************************************************************
set MODE  = $1
set SET   = $2
set DIR   = '~/Store/dstrtaunu/mcprod/gsim/'$MODE'/mdst_set'$SET'/'
set MDL    = 'DSTRTAUNU'
set TAG   = '_set'$SET    # _NAME
set FLAG  = 'sigMC'
set list = `awk '{print $1}' ~/script/nBB.txt`
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ${TAG} ]"

foreach f($list) 
  ./mk_script_sigMC.sh $DIR $MDL/$TAG $f $FLAG $MODE
end
