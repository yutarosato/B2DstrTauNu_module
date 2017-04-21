#! /bin/tcsh -f

if( $#argv != 2 )then
    echo " Usage : $0 [0p0,0p1,...1p0] [stream]"
    exit 1
endif

#****************************************************************************************************************
set MODE   = "Dstaunu"
set TYPE   = "Mixed"
set PAR    = $1
set STREAM = $2
set DIR   = "~/Store/dstrtaunu/mcprod_np_nagoya/store/BSemiTauonic_BSTD_2HDMTYPE2-${PAR}_${TYPE}_${MODE}/mdst/${STREAM}/"
set MDL    = 'DSTRTAUNU'
set FLAG  = 'npMC'
set list = `awk '{print $1}' ~/script/nBB.txt`
#****************************************************************************************************************
#echo "[ ${MDL} ${FLAG} ${PAR} ${STREAM}]"

foreach f($list) 
  ./mk_script_npMC.sh $DIR $MDL/$PAR/$STREAM $f $FLAG
end
