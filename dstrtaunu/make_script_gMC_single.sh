#! /bin/tcsh -f
if( $#argv != 1 )then
    echo " Usage : $0 [EVENT]"
    echo '[EVENT ] = {uds, charm, mixed, charged}'
    exit 1
endif

#****************************************************************************************************************
set MDL    = 'DSTRTAUNU'
set FLAG   = 'gMC'
set EVENT  = $1 # [uds, charm, mixed, charmed]
set STREAM = '0'
set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15/'
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set CNT=0
set POST=0
		   
foreach f($DIR"/"$EVENT"/gMC_"*".mdst")
     @ CNT += 1
     echo -n $CNT ' '
     ./mk_script_gMC_single.sh $MDL/$EVENT/$STREAM $f $FLAG $CNT
end
