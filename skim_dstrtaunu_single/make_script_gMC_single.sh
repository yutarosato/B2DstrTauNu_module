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
set list_data = `awk '{print $2}' ~syutaro/script/data_list_file_on_resonance_${EVENT}_e055_s0${STREAM}.txt`
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set CNT=0
set POST=0
		   
foreach f($list_data)
     @ CNT += 1
     echo -n $CNT ' '
     ./mk_script_gMC_single.sh $MDL/$EVENT/$STREAM $f $FLAG $CNT
end
