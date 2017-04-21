#! /bin/tcsh -f

if( $#argv != 2 )then
    echo " Usage : $0 [0p0/0p3/0p5/1p0] [set-name]"
    exit 1
endif

#****************************************************************************************************************
set MODE  = $1
set SET   = $2
set DIR   = '~/Store2/dstrtaunu/mcprod_np/evtgen/mdst_double/'$MODE'/' # sampels for double semileptonic decay events only for B0 -> D* tau nu
set MDL    = 'DSTRTAUNU'
set TAG   = $SET
set FLAG  = 'NPMC'
set list = `awk '{print $1}' ~/script/nBB.txt`
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ${TAG} ]"

foreach f($list) 
  ./mk_script_2HDM_double.sh $DIR $MDL/$TAG $f $FLAG $MODE
end
