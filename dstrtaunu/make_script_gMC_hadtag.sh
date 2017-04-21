#! /bin/tcsh -f

if( $#argv != 2 )then
    echo " Usage : $0 [EVENT][STREAM]"
    echo '[EVENT ] = {uds, charm, mixed, charged}'
    echo '[STREAM] = {0-5(9)}'
    exit 1
endif

#****************************************************************************************************************
set DIR    = '~/Store/dstrtaunu/modules/skim_dstrtaunu_hadtag/mdst_hadtag/mdst16/'
set MDL    = 'DSTRTAUNU'
set FLAG   = 'gMC'
set EVENT  = $1 # [uds, charm, mixed, charmed]
set STREAM = $2 # [ 0-5 ]
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set CNT=0

foreach f($DIR"/"$EVENT"/gMC_"$EVENT"_e"*s0$STREAM*".mdst")
     @ CNT += 1
     echo -n $CNT ' '
     ./mk_script_gMC_hadtag.sh $MDL/$EVENT/$STREAM $f $FLAG $CNT
end
