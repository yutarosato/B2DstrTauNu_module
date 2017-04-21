#! /bin/tcsh -f
#****************************************************************************************************************
#set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst6/'
#set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst14/'
set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst20/'
#set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst23/'

set MDL    = 'DSTRTAUNU'
set FLAG   = 'RD'
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set CNT=0

foreach f($DIR"/rd/RD_e"*".mdst")
     @ CNT += 1
     echo -n $CNT ' '
     ./mk_script_RD.sh $MDL $f $FLAG $CNT
end
