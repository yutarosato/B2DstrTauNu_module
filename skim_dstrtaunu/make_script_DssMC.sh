#! /bin/tcsh -f

#****************************************************************************************************************
set list   = '/home/belle/syutaro/script/file_list_DssMC.txt'
set MDL    = 'DSTRTAUNU'
set FLAG  = 'DssMC'
#****************************************************************************************************************
echo "[ ${MDL} ${FLAG} ]"

set NLINE = `wc -l ${list} | cut -f1 -d " "`
set LIMIT = 16
set CNT_START = 1
set CNT_END = 1
set TOTAL_CNT = 1
set CNT = 1  # for count limit
set CNT2 = 1 # for c00x
set postexpNo = ""
foreach f(`cat $list`)
  set expNo  = `basename $f | cut -d "_" -f 3`

   #echo -n "${expNo} ${postexpNo} ${CNT} ${TOTAL_CNT} ${f} " 

  if( ( (${expNo} != ${postexpNo}) || (${CNT} == ${LIMIT}) || ${TOTAL_CNT} == ${NLINE} ) && ${TOTAL_CNT} != '1'  ) then
    @ CNT = 1
	@ CNT_END = $TOTAL_CNT - 1
      if( ${TOTAL_CNT} == ${NLINE} ) then
	@ CNT_END += 1
      endif 
    ./mk_script_DssMC.sh $MDL $FLAG $postexpNo $CNT2 $CNT_START $CNT_END $list
    #echo -n "${CNT_START} ${CNT_END}"

    @ CNT2 += 1
    if( ${expNo} != ${postexpNo} ) then
	@ CNT2 = 1
    endif
    @ CNT_START = $CNT_END + 1
  endif

   #echo ""
   set postexpNo = $expNo
   @ CNT  += 1
   @ TOTAL_CNT += 1
end
