#! /bin/tcsh -f

#****************************************************************************************************************
set MDL    = 'DSTRTAUNU'
set FLAG   = 'RD'
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
     ./mk_script_RD.sh $MDL $f $FLAG $CNT
end
