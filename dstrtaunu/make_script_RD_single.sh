#! /bin/tcsh -f

#****************************************************************************************************************
set MDL    = 'DSTRTAUNU'
set FLAG   = 'RD'
set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15/rd'
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set CNT=0
set POST=0
		   
foreach f($DIR"/RD_"*".mdst")
     @ CNT += 1
     echo -n $CNT ' '
     ./mk_script_RD_single.sh $MDL $f $FLAG $CNT
end
