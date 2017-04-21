#! /bin/tcsh -f

set DIR    = $1
set MDL    = `echo $2 | cut -d "/" -f 1`
set TAG    = `echo $2 | cut -d "/" -f 2`
set EXP    = $3
set FLAG   = $4
set MODE   = $5
set TYPE   = $6
      
set  OUT = "script/OP_${TYPE}_${MODE}_double/${FLAG}_${TYPE}_${MODE}_${EXP}_set${TAG}.sh"
set  LOG =    "log/OP_${TYPE}_${MODE}_double/${FLAG}_${TYPE}_${MODE}_${EXP}_set${TAG}.log"
set  HBK =    "hbk/OP_${TYPE}_${MODE}_double/${FLAG}_${TYPE}_${MODE}_${EXP}_set${TAG}.hbk"
     
mkdir -p "script/OP_${TYPE}_${MODE}_double"
#****************************************
      
echo '#! /bin/tcsh -f\n' > $OUT

echo 'setenv MY_TOP_DIR $HOME/dstrtaunu/modules/np_dstrtaunu_gen' >> $OUT
echo 'setenv BELLE_MESSAGE_LEVEL INFO' >> $OUT
echo 'setenv USE_GRAND_REPROCESS_DATA 1' >> $OUT
echo 'setenv FPDA_TIMEOUT 0' >> $OUT
echo 'source /sw/belle/local/etc/cshrc_general\n' >> $OUT # KEKCC
#echo 'source /belle/local/etc/cshrc_general\n' >> $OUT # NAGOYA


echo 'echo    $LD_LIBRARY_PATH'     >> $OUT
echo 'echo    $BASF_MODULE_DIR'     >> $OUT
echo 'echo    $PANTHER_TABLE_DIR\n' >> $OUT

echo "mkdir -p    log/OP_${TYPE}_${MODE}_double" >> $OUT
echo "mkdir -p    hbk/OP_${TYPE}_${MODE}_double" >> $OUT

echo 'basf << EOF > ' $LOG'\n' >> $OUT

echo 'nprocess set 0\n' >> $OUT
  
echo 'path create main'     >> $OUT
echo 'path create analysis\n' >> $OUT

echo 'path add_module main' $MDL           >> $OUT
echo 'path add_condition main >:0:EXIT'    >> $OUT
echo 'path add_condition main <=:0:KILL\n' >> $OUT

echo '' >> $OUT

echo 'initialize' >> $OUT
echo 'histogram define '$HBK'\n' >> $OUT

set CNT = 0
foreach f ($DIR"evtgen_exp_"$EXP"_"*"__"$TAG"_"*.mdst)
  echo 'process_event '$f >> $OUT
  @ CNT++
end

echo '\nterminate' >> $OUT
echo 'EOF\n' >> $OUT

echo 'h2root '$HBK '&& rm '$HBK >> $OUT 

echo $CNT $OUT

chmod 755 $OUT
