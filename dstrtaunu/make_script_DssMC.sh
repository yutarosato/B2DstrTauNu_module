#! /bin/tcsh -f

#****************************************************************************************************************
#set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst20/'
set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst23/'
set MDL    = 'DSTRTAUNU'
set FLAG  = 'DssMC'
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set CNT=0

foreach f($DIR"/DssMC/"*".mdst")
     @ CNT += 1
     ./mk_script_DssMC.sh $MDL $f $FLAG $CNT
end
