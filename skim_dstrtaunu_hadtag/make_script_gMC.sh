#! /bin/tcsh -f

if( $#argv != 2 )then
    echo " Usage : $0 [EVENT][STREAM]"
    echo '[EVENT ] = {uds, charm, mixed, charged}'
    echo '[STREAM] = {0-5(9)}'
    exit 1
endif

#****************************************************************************************************************
set MDL    = 'DSTRTAUNU'
set FLAG   = 'gMC'
set EVENT  = $1 # [uds, charm, mixed, charmed]
set STREAM = $2 # [ 0-5 ]
set list_data = `awk '{print $1}'  ~syutaro/script/data_list_on_resonance.txt`  
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set CNT=0
set POST=0
		   
foreach f($list_data)
     set expNo  = `echo $f | cut -d "/" -f 1`
     if( $expNo != $POST ) then
       @ CNT = 0
     endif
     @ POST = $expNo
     @ CNT += 1
     echo -n $CNT ' '
     ./mk_script_gMC.sh $MDL/$EVENT/$STREAM $f $FLAG $CNT
end
