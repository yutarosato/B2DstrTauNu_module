#! /bin/tcsh -f

# 26 expNo, 21177 mdst files
#****************************************************************************************************************
set MDL    = 'DSTRTAUNU'
set FLAG   = 'RD'
#****************************************************************************************************************

foreach No (`cat ~syutaro/script/nBB2.txt`)
set expNo = `echo $No | sed -e "s|^0||"`
echo -n $expNo" "
  ~syutaro/script/find_rd $expNo 1 9999 on_resonance | awk '{print $2}' > tmp.list.rd
  set CNT=0
  foreach f(`cat tmp.list.rd`)
    @ CNT += 1
    ./mk_script_RD_onebyone.sh $MDL $f $expNo $FLAG $CNT
  end
  echo $CNT
end

