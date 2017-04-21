#! /bin/tcsh -f

#****************************************************************************************************************
set MDL    = 'DSTRTAUNU'
set FLAG   = 'RD'
set list_data = `awk '{print $2}' ~syutaro/script/data_list_file_on_resonance_rd_e055.txt`
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set CNT=0
set POST=0
		   
foreach f($list_data)
     @ CNT += 1
     echo -n $CNT ' '
     ./mk_script_RD_single.sh $MDL $f $FLAG $CNT
end
