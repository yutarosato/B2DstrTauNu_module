#! /bin/tcsh -f

if( $#argv != 3 )then
    echo " Usage : $0  1p05 S1 [set-name]"
    echo " Usage : $0 m2p10 V1 [set-name]"
    echo " Usage : $0  1p88 V2 [set-name]"
    exit 1
endif

#****************************************************************************************************************
set MODE  = $1
set TYPE  = $2
set SET   = $3
set DIR   = "~/Store2/dstrtaunu/mcprod_np/evtgen_npoperator/${TYPE}/mdst/${MODE}/"  # samples for all modes
set MDL    = 'DSTRTAUNU'
set TAG   = $SET
set FLAG  = 'NPMC'
set list = `awk '{print $1}' ~/script/nBB.txt`
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ${TAG} ${MODE} ${TYPE} ]"

foreach f($list) 
  ./mk_script_NPMC.sh $DIR $MDL/$TAG $f $FLAG $MODE $TYPE
end
