#! /bin/tcsh -f

#****************************************************************************************************************
set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu_hadtag/mdst_hadtag/mdst16/'
set MDL    = 'DSTRTAUNU'
set FLAG   = 'RD'
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set CNT=0

foreach f($DIR"/rd/RD_e"*".mdst")
     @ CNT += 1
     echo -n $CNT ' '
     ./mk_script_RD_hadtag.sh $MDL $f $FLAG $CNT
end
